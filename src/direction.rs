//! This module contains the Direction struct, which represents a normalized vector in 3D space.

use serde::ser::SerializeTuple;
use serde::Serializer;
use serde::{Deserialize, Serialize};
use simple_si_units::{base::Distance, geometry::Angle};
use std::fmt::Display;
use std::ops::Neg;

use crate::equatorial::Equatorial;
use crate::error::AstroCoordsError;
use crate::traits::*;
use crate::transformations::rotations::*;
use crate::{angle_helper::*, NORMALIZATION_THRESHOLD};

use super::{
    cartesian::Cartesian, earth_equatorial::EarthEquatorial, ecliptic::Ecliptic,
    spherical::Spherical,
};

/// The Direction struct represents a normalised vector in 3D space.
#[derive(Debug, Clone, PartialEq)]
pub struct Direction {
    pub(super) x: f64,
    pub(super) y: f64,
    pub(super) z: f64,
}

impl Mathematical for Direction {}

impl Direction {
    /// A normalised vector pointing in x-Direction.
    pub const X: Direction = Direction {
        x: 1.,
        y: 0.,
        z: 0.,
    };

    /// A normalised vector pointing in y-Direction.
    pub const Y: Direction = Direction {
        x: 0.,
        y: 1.,
        z: 0.,
    };

    /// A normalised vector pointing in z-Direction.
    pub const Z: Direction = Direction {
        x: 0.,
        y: 0.,
        z: 1.,
    };

    const SERIALIZATION_ACCURACY: f64 = 1e-3;

    /// Creates a new Direction from the given coordinates.
    ///
    /// If the length of the vector is below the `NORMALIZATION_THRESHOLD`, an error is returned.
    pub fn new(x: f64, y: f64, z: f64) -> Result<Self, AstroCoordsError> {
        let length = (x * x + y * y + z * z).sqrt();
        if length < NORMALIZATION_THRESHOLD {
            Err(AstroCoordsError::NormalizingZeroVector)
        } else {
            Ok(Direction {
                x: x / length,
                y: y / length,
                z: z / length,
            })
        }
    }

    /// Returns the Direction as an array.
    pub fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Creates a Direction from an array.
    pub fn from_array(array: [f64; 3]) -> Result<Self, AstroCoordsError> {
        Direction::new(array[0], array[1], array[2])
    }

    /// Returns a Cartesian struct with the specified length, pointing in the direction.
    ///
    /// # Example
    /// ```
    /// use simple_si_units::base::Distance;
    /// use astro_coords::{direction::Direction, cartesian::Cartesian};
    ///
    /// let direction = Direction::new(1., 1., 1.).unwrap();
    /// let length = Distance::from_meters(10.);
    /// let ordinate_length = length / 3f64.sqrt();
    /// let cartesian = direction.to_cartesian(length);
    /// let expected = Cartesian::new(ordinate_length, ordinate_length, ordinate_length);
    /// assert!(cartesian.eq_within(&expected, Distance::from_meters(1e-5)));
    /// ```
    pub fn to_cartesian(&self, length: Distance<f64>) -> Cartesian {
        Cartesian::new(self.x * length, self.y * length, self.z * length)
    }

    /// Returns the spherical coordinates of the Direction.
    ///
    /// # Example
    /// ```
    /// use astro_coords::direction::Direction;
    ///
    /// let direction = Direction::X;
    /// let spherical = direction.to_spherical();
    /// assert!((spherical.latitude.to_degrees() - 0.).abs() < 1e-5);
    /// assert!((spherical.longitude.to_degrees() - 0.).abs() < 1e-5);
    ///
    /// let direction = Direction::Y;
    /// let spherical = direction.to_spherical();
    /// assert!((spherical.latitude.to_degrees() - 0.).abs() < 1e-5);
    /// assert!((spherical.longitude.to_degrees() - 90.).abs() < 1e-5);
    ///
    /// let direction = Direction::Z;
    /// let spherical = direction.to_spherical();
    /// assert!((spherical.latitude.to_degrees() - 90.).abs() < 1e-5);
    /// assert!((spherical.longitude.to_degrees() - 0.).abs() < 1e-5);
    /// ```
    pub fn to_spherical(&self) -> Spherical {
        // Direction is normalised and thus guaranteed to produce valid spherical coordinates.
        Spherical::cartesian_to_spherical((self.x, self.y, self.z))
            .unwrap_or(Spherical::new(ANGLE_ZERO, ANGLE_ZERO))
    }

    /// Returns the x-ordinate of the Direction.
    pub fn x(&self) -> f64 {
        self.x
    }

    /// Returns the y-ordinate of the Direction.
    pub fn y(&self) -> f64 {
        self.y
    }

    /// Returns the z-ordinate of the Direction.
    pub fn z(&self) -> f64 {
        self.z
    }

    /// Returns true if the Direction is equal to the other Direction within the specified accuracy.
    ///
    /// # Example
    /// ```
    /// use astro_coords::direction::Direction;
    ///
    /// let direction = Direction::new(1., 1., 1.).unwrap();
    /// let other = Direction::new(1., 1., 1.).unwrap();
    /// assert!(direction.eq_within(&other, 1e-5));
    /// ```
    pub fn eq_within(&self, other: &Direction, accuracy: f64) -> bool {
        (self.x - other.x).abs() < accuracy
            && (self.y - other.y).abs() < accuracy
            && (self.z - other.z).abs() < accuracy
    }

    /// Returns the angle between the Direction and the other Direction.
    ///
    /// # Example
    /// ```
    /// use astro_coords::direction::Direction;
    ///
    /// let direction = Direction::X;
    /// let other = Direction::Y;
    /// let angle = direction.angle_to(&other);
    /// assert!((angle.to_degrees() - 90.).abs() < 1e-5);
    /// ```
    pub fn angle_to(&self, other: &Direction) -> Angle<f64> {
        let (ax, ay, az) = (self.x(), self.y(), self.z());
        let (bx, by, bz) = (other.x(), other.y(), other.z());

        let cosine_argument = ax * bx + ay * by + az * bz; //Directions have unit Distance
        safe_acos(cosine_argument)
    }

    /// Constructs a Direction that is guaranteed to be orthogonal to the Direction.
    ///
    /// # Example
    /// ```
    /// use astro_coords::direction::Direction;
    ///
    /// let direction = Direction::new(1., 2., 3.).unwrap();
    /// let orthogonal = direction.some_orthogonal_vector();
    /// assert!(orthogonal.angle_to(&direction).to_degrees() - 90. < 1e-5);
    /// ```
    pub fn some_orthogonal_vector(&self) -> Direction {
        let ortho = if self.x().abs() > NORMALIZATION_THRESHOLD {
            self.cross_product(&Self::Y)
        } else if self.y().abs() > NORMALIZATION_THRESHOLD {
            self.cross_product(&Self::Z)
        } else {
            self.cross_product(&Self::X)
        };
        ortho.unwrap_or(Self::Z)
    }

    /// Constructs the cross product of the Direction and the other Direction.
    ///
    /// # Example
    /// ```
    /// use astro_coords::direction::Direction;
    ///
    /// let direction = Direction::X;
    /// let other = Direction::Y;
    /// let cross = direction.cross_product(&other).unwrap();
    /// assert!(cross.eq_within(&Direction::Z, 1e-5));
    /// ```
    pub fn cross_product(&self, other: &Direction) -> Result<Direction, AstroCoordsError> {
        let (ax, ay, az) = (self.x, self.y, self.z);
        let (bx, by, bz) = (other.x(), other.y(), other.z());

        let cx = ay * bz - az * by;
        let cy = az * bx - ax * bz;
        let cz = ax * by - ay * bx;

        Direction::new(cx, cy, cz)
    }

    /// Calculates the dot product of the Direction and the other Direction.
    ///
    /// # Example
    /// ```
    /// use astro_coords::direction::Direction;
    ///
    /// let direction = Direction::X;
    /// let other = Direction::Y;
    /// let dot_product = direction.dot_product(&other);
    /// assert!((dot_product - 0.).abs() < 1e-5);
    /// ```
    pub fn dot_product(&self, other: &Direction) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn to_earth_equatorial(&self) -> EarthEquatorial {
        let dir_in_equatorial = self.rotated(EARTH_AXIS_TILT, &Direction::X);
        let spherical = dir_in_equatorial.to_spherical();
        EarthEquatorial::new(spherical.longitude, spherical.latitude)
    }

    pub fn to_ecliptic(&self) -> Ecliptic {
        Ecliptic::new(self.to_spherical())
    }

    pub fn to_equatorial(&self, axis: Direction) -> Equatorial {
        let dir_in_equatorial = self.passive_rotation_to_new_z_axis(&axis);
        let spherical = dir_in_equatorial.to_spherical();
        Equatorial::new(spherical, axis)
    }
}

impl ActiveRotation<Direction> for Direction {
    /// Returns a new Direction that is rotated by the specified angle around the specified axis.
    ///
    /// # Example
    /// ```
    /// use simple_si_units::geometry::Angle;
    /// use astro_coords::{direction::Direction, traits::*};
    ///
    /// let direction = Direction::X;
    /// let angle = Angle::from_degrees(90.);
    /// let axis = Direction::Z;
    /// let rotated = direction.rotated(angle, &axis);
    /// assert!(rotated.eq_within(&Direction::Y, 1e-5));
    /// ```
    fn rotated(&self, angle: Angle<f64>, axis: &Direction) -> Direction {
        let (x, y, z) = rotated_tuple((self.x, self.y, self.z), angle, axis);
        Direction { x, y, z }
    }

    /// Returns a new Direction that is rotated around the x-axis by a certain angle.
    ///
    /// This is an ever so tiny bit faster than `rotated()`.
    ///
    /// # Example
    /// ```
    /// use simple_si_units::geometry::Angle;
    /// use astro_coords::{direction::Direction, traits::*};
    ///
    /// let direction = Direction::Y;
    /// let angle = Angle::from_degrees(90.);
    /// let rotated = direction.rotated_x(angle);
    /// assert!(rotated.eq_within(&Direction::Z, 1e-5));
    /// ```
    fn rotated_x(&self, angle: Angle<f64>) -> Direction {
        let (x, y, z) = rotated_x_tuple((self.x, self.y, self.z), angle);
        Direction { x, y, z }
    }

    /// Returns a new Direction that is rotated around the y-axis by a certain angle.
    ///
    /// This is an ever so tiny bit faster than `rotated()`.
    ///
    /// # Example
    /// ```
    /// use simple_si_units::geometry::Angle;
    /// use astro_coords::{direction::Direction, traits::*};
    ///
    /// let direction = Direction::Z;
    /// let angle = Angle::from_degrees(90.);
    /// let rotated = direction.rotated_y(angle);
    /// assert!(rotated.eq_within(&Direction::X, 1e-5));
    /// ```
    fn rotated_y(&self, angle: Angle<f64>) -> Direction {
        let (x, y, z) = rotated_y_tuple((self.x, self.y, self.z), angle);
        Direction { x, y, z }
    }

    /// Returns a new Direction that is rotated around the z-axis by a certain angle.
    ///
    /// This is an ever so tiny bit faster than `rotated()`.
    ///
    /// # Example
    /// ```
    /// use simple_si_units::geometry::Angle;
    /// use astro_coords::{direction::Direction, traits::*};
    ///
    /// let direction = Direction::X;
    /// let angle = Angle::from_degrees(90.);
    /// let rotated = direction.rotated_z(angle);
    /// assert!(rotated.eq_within(&Direction::Y, 1e-5));
    /// ```
    fn rotated_z(&self, angle: Angle<f64>) -> Direction {
        let (x, y, z) = rotated_z_tuple((self.x, self.y, self.z), angle);
        Direction { x, y, z }
    }

    /// Returns the Direction that results from actively rotating the Direction to the new z-axis, in a manner that preserves the old z-projection of the x-axis.
    ///
    /// This method is for example used to convert from equatorial coordinates to ecliptic coordinates.
    /// It operates in the following way:
    /// 1. The vector is rotated around the old x-axis by the angle between new and old z-axis.
    /// 2. The vector is rotated around the old z-axis by the angle between the new z-axis and the old y-axis, projected onto the old x-y plane.
    ///
    /// This is the inverse operation of `passive_rotation_to_new_z_axis`. See there for a somewhat intuitive example.
    fn active_rotation_to_new_z_axis(&self, new_z: &Direction) -> Direction {
        let (angle_to_old_z, polar_rotation_angle) =
            get_angle_to_old_z_and_polar_rotation_angle(new_z);
        self.rotated(-angle_to_old_z, &Self::X)
            .rotated(-polar_rotation_angle, &Self::Z)
    }
}

impl PassiveRotation<Direction> for Direction {
    /// Returns the Direction that results from passively rotating the Direction to the new z-axis, in a manner that preserves the old z-projection of the x-axis.
    ///
    /// This method is for example used to convert from ecliptic coordinates to equatorial coordinates.
    /// It operates in the following way:
    /// 1. The vector is rotated around the old z-axis by the angle between the new z-axis and the old y-axis, projected onto the old x-y plane.
    /// 2. The vector is rotated around the old x-axis by the angle between new and old z-axis.
    ///
    /// This is the inverse operation of `active_rotation_to_new_z_axis`.
    ///
    /// # Example
    /// ```
    /// use simple_si_units::geometry::Angle;
    /// use astro_coords::{direction::Direction, traits::*};
    ///
    /// // Suppose it is summer solstice and the sun is in y-direction in the ecliptic coordinate system.
    /// let dir_of_sun_in_ecliptic = Direction::Y;
    ///
    /// // Now we want to express the sun's direction in earth equatorial coordinates.
    /// // The rotation axis of the earth expressed in ecliptic coordinates is given by:
    /// let earth_axis_tilt = Angle::from_degrees(23.44);
    /// let earth_rotation_axis_in_ecliptic = Direction::Z.rotated_x(-earth_axis_tilt);
    ///
    /// // The sun's direction in earth equatorial coordinates is then:
    /// let dir_of_sun_in_equatorial = dir_of_sun_in_ecliptic.passive_rotation_to_new_z_axis(&earth_rotation_axis_in_ecliptic);
    ///
    /// // At summer solstice, the sun is highest in the sky in the northern hemisphere, so its x-projection is zero, and its y- and z-projection are both positive.
    /// println!("{}", dir_of_sun_in_equatorial);
    /// assert!(dir_of_sun_in_equatorial.x().abs() < 1e-5);
    /// assert!(dir_of_sun_in_equatorial.y() > 0.);
    /// assert!(dir_of_sun_in_equatorial.z() > 0.);
    /// ```
    fn passive_rotation_to_new_z_axis(&self, new_z: &Direction) -> Direction {
        let (angle_to_old_z, polar_rotation_angle) =
            get_angle_to_old_z_and_polar_rotation_angle(new_z);
        self.rotated(polar_rotation_angle, &Self::Z)
            .rotated(angle_to_old_z, &Self::X)
    }
}

fn get_angle_to_old_z_and_polar_rotation_angle(new_z: &Direction) -> (Angle<f64>, Angle<f64>) {
    let angle_to_old_z = new_z.angle_to(&Direction::Z);

    let axis_projected_onto_xy_plane = Direction::new(new_z.x(), new_z.y(), 0.);
    let mut polar_rotation_angle = ANGLE_ZERO;
    if let Ok(axis_projected_onto_xy_plane) = axis_projected_onto_xy_plane {
        polar_rotation_angle = axis_projected_onto_xy_plane.angle_to(&Direction::Y);
        if axis_projected_onto_xy_plane.x() < 0. {
            polar_rotation_angle = -polar_rotation_angle;
        }
    }
    (angle_to_old_z, polar_rotation_angle)
}

impl Neg for &Direction {
    type Output = Direction;

    fn neg(self) -> Self::Output {
        Direction {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Display for Direction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({:.2}, {:.2}, {:.2})", self.x, self.y, self.z)
    }
}

impl Serialize for Direction {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let array = self.to_array();
        let mut tuple_serializer = serializer.serialize_tuple(3)?;
        for value in &array {
            let value =
                (value / Self::SERIALIZATION_ACCURACY).round() * Self::SERIALIZATION_ACCURACY;
            tuple_serializer.serialize_element(&value)?;
        }
        tuple_serializer.end()
    }
}

impl<'de> Deserialize<'de> for Direction {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let array = <[f64; 3]>::deserialize(deserializer)?;
        Direction::from_array(array).map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_ACCURACY: f64 = 1e-5;
    const ROTATION_ACCURACY: Angle<f64> = Angle { rad: 1e-3 }; //Accos is a bit unstable

    #[test]
    fn from_spherical() {
        let values = vec![-1., 1., 10.];
        for x in values.iter() {
            for y in values.iter() {
                for z in values.iter() {
                    let x = Distance::from_meters(*x);
                    let y = Distance::from_meters(*y);
                    let z = Distance::from_meters(*z);
                    let cartesian = Cartesian::new(x, y, z);
                    let length = cartesian.length();
                    let expected_x = x / length;
                    let expected_y = y / length;
                    let expected_z = z / length;
                    let direction = cartesian.to_spherical().unwrap().to_direction();

                    assert!((direction.x() - expected_x).abs() < TEST_ACCURACY);
                    assert!((direction.y() - expected_y).abs() < TEST_ACCURACY);
                    assert!((direction.z() - expected_z).abs() < TEST_ACCURACY);
                }
            }
        }
    }

    #[test]
    fn angle_between_is_half_turn() {
        let expected: Angle<f64> = HALF_CIRC;

        let angle = Direction::X.angle_to(&(-&Direction::X));
        println!("angle: {}", angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));

        let angle = Direction::Y.angle_to(&(-&Direction::Y));
        println!("angle: {}", angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));

        let angle = Direction::Z.angle_to(&(-&Direction::Z));
        println!("angle: {}", angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));

        let angle1 = Direction::new(1., 1., 0.).unwrap();
        let angle2 = Direction::new(-1., -1., 0.).unwrap();
        let angle = angle1.angle_to(&angle2);
        println!("angle: {}", angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));

        let angle1 = Direction::new(1., 0., 1.).unwrap();
        let angle2 = Direction::new(-1., 0., -1.).unwrap();
        let angle = angle1.angle_to(&angle2);
        println!("angle: {}", angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));

        let angle1 = Direction::new(0., 1., 1.).unwrap();
        let angle2 = Direction::new(0., -1., -1.).unwrap();
        let angle = angle1.angle_to(&angle2);
        println!("angle: {}", angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));
    }

    #[test]
    fn angle_between_is_quarter_turn() {
        let expected = QUARTER_CIRC;

        let angle = Direction::X.angle_to(&Direction::Y);
        println!("expected: {}, actual: {}", expected, angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));

        let angle = Direction::X.angle_to(&Direction::Z);
        println!("expected: {}, actual: {}", expected, angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));

        let angle = Direction::Y.angle_to(&Direction::Z);
        println!("expected: {}, actual: {}", expected, angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));

        let angle1 = Direction::new(1., 1., 0.).unwrap();
        let angle2 = Direction::new(1., -1., 0.).unwrap();
        let angle = angle1.angle_to(&angle2);
        println!("expected: {}, actual: {}", expected, angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));

        let angle1 = Direction::new(1., 0., 1.).unwrap();
        let angle2 = Direction::new(1., 0., -1.).unwrap();
        let angle = angle1.angle_to(&angle2);
        println!("expected: {}, actual: {}", expected, angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));
    }

    #[test]
    fn angle_between_is_zero() {
        let expected: Angle<f64> = Angle::from_rad(0.);

        let angle = Direction::X.angle_to(&Direction::X);
        println!("expected: {}, actual: {}", expected, angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));

        let angle = Direction::Y.angle_to(&Direction::Y);
        println!("expected: {}, actual: {}", expected, angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));

        let angle = Direction::Z.angle_to(&Direction::Z);
        println!("expected: {}, actual: {}", expected, angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));

        let angle1 = Direction::new(1., 1., 0.).unwrap();
        let angle2 = Direction::new(1., 1., 0.).unwrap();
        let angle = angle1.angle_to(&angle2);
        println!("expected: {}, actual: {}", expected, angle);
        assert!(angle_eq_within(angle, expected, ROTATION_ACCURACY));
    }

    #[test]
    fn test_cross_product() {
        let angle1 = Direction::new(1., 0., 0.).unwrap();
        let angle2 = Direction::new(0., 1., 0.).unwrap();
        let expected = Direction::new(0., 0., 1.).unwrap();
        let actual = angle1.cross_product(&angle2).unwrap();
        println!("expected: {}, actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));

        let angle1 = Direction::new(1., 0., 0.).unwrap();
        let angle2 = Direction::new(0., 0., 1.).unwrap();
        let expected = Direction::new(0., -1., 0.).unwrap();
        let actual = angle1.cross_product(&angle2).unwrap();
        println!("expected: {}, actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));

        let angle1 = Direction::new(0., 1., 0.).unwrap();
        let angle2 = Direction::new(0., 0., 1.).unwrap();
        let expected = Direction::new(1., 0., 0.).unwrap();
        let actual = angle1.cross_product(&angle2).unwrap();
        println!("expected: {}, actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));

        let angle1 = Direction::new(0., 1., 0.).unwrap();
        let angle2 = Direction::new(0., 0., -1.).unwrap();
        let expected = Direction::new(-1., 0., 0.).unwrap();
        let actual = angle1.cross_product(&angle2).unwrap();
        println!("expected: {}, actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));

        let angle1 = Direction::new(0., 0., 1.).unwrap();
        let angle2 = Direction::new(0., 1., 0.).unwrap();
        let expected = Direction::new(-1., 0., 0.).unwrap();
        let actual = angle1.cross_product(&angle2).unwrap();
        println!("expected: {}, actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));

        let angle1 = Direction::new(0., 0., 1.).unwrap();
        let angle2 = Direction::new(1., 0., 0.).unwrap();
        let expected = Direction::new(0., 1., 0.).unwrap();
        let actual = angle1.cross_product(&angle2).unwrap();
        println!("expected: {}, actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));
    }

    #[test]
    fn cross_product_is_always_orthogonal() {
        let ordinates = vec![-1., 0., 1., 10.];
        for x in ordinates.clone().iter() {
            for y in ordinates.clone().iter() {
                for z in ordinates.clone().iter() {
                    for u in ordinates.clone().iter() {
                        for v in ordinates.clone().iter() {
                            for w in ordinates.clone().iter() {
                                let a = Direction::new(*x, *y, *z);
                                let b = Direction::new(*u, *v, *w);
                                if a.is_err() || b.is_err() {
                                    continue;
                                }
                                let a = a.unwrap();
                                let b = b.unwrap();
                                println!("a: {}, b: {}", a, b);
                                let cross = a.cross_product(&b);
                                if a.eq_within(&b, TEST_ACCURACY)
                                    || a.eq_within(&-&b, TEST_ACCURACY)
                                {
                                    assert!(cross.is_err());
                                } else {
                                    let cross = cross.unwrap();
                                    let overlap_with_a = cross.dot_product(&a);
                                    let overlap_with_b = cross.dot_product(&b);
                                    assert!(overlap_with_a.abs() < TEST_ACCURACY);
                                    assert!(overlap_with_b.abs() < TEST_ACCURACY);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_revertability_of_rotation() {
        let ordinates: Vec<f64> = vec![-1., 0., 1., 10.];
        for x1 in ordinates.clone() {
            for y1 in ordinates.clone() {
                for z1 in ordinates.clone() {
                    for x2 in ordinates.clone() {
                        for y2 in ordinates.clone() {
                            for z2 in ordinates.clone() {
                                let original_dir = Direction::new(x1, y1, z1);
                                let new_z_axis = Direction::new(x2, y2, z2);
                                if original_dir.is_err() || new_z_axis.is_err() {
                                    continue;
                                }
                                let original_dir = original_dir.unwrap();
                                let new_z_axis = new_z_axis.unwrap();

                                let after_active_rotation =
                                    original_dir.active_rotation_to_new_z_axis(&new_z_axis);
                                let transformed_back = after_active_rotation
                                    .passive_rotation_to_new_z_axis(&new_z_axis);

                                println!("new_z_axis: {}", new_z_axis);
                                println!("original_dir: {}", original_dir);
                                println!("after_active_rotation: {}", after_active_rotation);
                                println!("transformed_back: {}", transformed_back);
                                assert!(original_dir.eq_within(&transformed_back, TEST_ACCURACY));
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_north_after_active_rotation() {
        let ordinates: Vec<f64> = vec![-1., 0., 1., 10.];
        for x in ordinates.clone() {
            for y in ordinates.clone() {
                for z in ordinates.clone() {
                    let new_z_axis = Direction::new(x, y, z);
                    if new_z_axis.is_err() {
                        continue;
                    }
                    let new_z_axis = new_z_axis.unwrap();

                    let expected = new_z_axis.clone();
                    let actual = Direction::Z.active_rotation_to_new_z_axis(&new_z_axis);

                    println!("new_z_axis: {}", new_z_axis);
                    println!("expected: {}", expected);
                    println!("actual: {}", actual);
                    assert!(expected.eq_within(&actual, TEST_ACCURACY));
                }
            }
        }
    }

    #[test]
    fn active_rotation_to_z_changes_nothing() {
        let ordinates: Vec<f64> = vec![-1., 0., 1., 10.];
        for x in ordinates.clone() {
            for y in ordinates.clone() {
                for z in ordinates.clone() {
                    let dir = Direction::new(x, y, z);
                    if dir.is_err() {
                        continue;
                    }
                    let dir = dir.unwrap();

                    let expected = dir.clone();
                    let actual = dir.active_rotation_to_new_z_axis(&Direction::Z);

                    println!("expected: {}", expected);
                    println!("actual: {}", actual);
                    assert!(expected.eq_within(&actual, TEST_ACCURACY));
                }
            }
        }
    }

    #[test]
    fn test_north_after_passive_rotation() {
        let ordinates: Vec<f64> = vec![-1., 0., 1., 10.];
        for x in ordinates.clone() {
            for y in ordinates.clone() {
                for z in ordinates.clone() {
                    let new_z_axis = Direction::new(x, y, z);
                    if new_z_axis.is_err() {
                        continue;
                    }
                    let new_z_axis = new_z_axis.unwrap();

                    let expected = Direction::Z;
                    let actual = new_z_axis.passive_rotation_to_new_z_axis(&new_z_axis);

                    println!("new_z_axis: {}", new_z_axis);
                    println!("expected: {}", expected);
                    println!("actual: {}", actual);
                    assert!(expected.eq_within(&actual, TEST_ACCURACY));
                }
            }
        }
    }

    #[test]
    fn direction_from_large_vector_is_ok() {
        let x = Distance::from_lyr(2000.);
        let y = Distance::from_lyr(1e-10);
        let z = Distance::from_lyr(-2000.);
        let cartesian = Cartesian::new(x, y, z);
        let expected = Direction::new(1., 0., -1.).unwrap();
        let actual = cartesian.to_direction().unwrap();
        println!("expected: {}, actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));
    }

    #[test]
    fn angle_between_two_close_directions() {
        // The maths is correct, but this is really unstable due to accos.
        let test_accuracy = Angle { rad: 0.03 };
        let a = Direction {
            x: -0.085366555,
            y: -0.8412673,
            z: -0.5338369,
        };
        let b = Direction {
            x: -0.08536589,
            y: -0.84126735,
            z: -0.5338369,
        };
        let angle = a.angle_to(&b);
        println!("angle: {}", angle);
        assert!(angle_eq_within(angle, ANGLE_ZERO, test_accuracy));
    }

    #[test]
    fn serialization() {
        let dir = Direction::new(1.23, -0.01, 1e-8).unwrap();
        let serialized = serde_json::to_string(&dir).unwrap();
        println!("{:?}", dir);
        println!("{}", serialized);
        assert_eq!(serialized, "[1.0,-0.008,0.0]");

        let deserialized: Direction = serde_json::from_str(&serialized).unwrap();
        println!("{:?}", deserialized);
        assert!(deserialized.eq_within(&dir, Direction::SERIALIZATION_ACCURACY));
    }
}
