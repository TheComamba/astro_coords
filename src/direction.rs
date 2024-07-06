use serde::ser::SerializeTuple;
use serde::Serializer;
use serde::{Deserialize, Serialize};
use simple_si_units::{base::Distance, geometry::Angle};
use std::fmt::Display;
use std::ops::Neg;

use crate::angle_helper::*;
use crate::equatorial::EquatorialCoordinates;
use crate::error::AstroCoordsError;

use super::{
    cartesian::CartesianCoordinates, earth_equatorial::EarthEquatorialCoordinates,
    ecliptic::EclipticCoordinates, spherical::SphericalCoordinates,
    transformations::rotations::rotated_tuple,
};

pub(super) const NORMALIZATION_THRESHOLD: f64 = 1e-5;

#[derive(Debug, Clone, PartialEq)]
pub struct Direction {
    pub(super) x: f64,
    pub(super) y: f64,
    pub(super) z: f64,
}

impl Direction {
    const SERIALIZATION_ACCURACY: f64 = 1e-3;

    pub fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    pub fn from_array(array: [f64; 3]) -> Result<Self, AstroCoordsError> {
        Direction::new(array[0], array[1], array[2])
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

impl Direction {
    pub const X: Direction = Direction {
        x: 1.,
        y: 0.,
        z: 0.,
    };
    pub const Y: Direction = Direction {
        x: 0.,
        y: 1.,
        z: 0.,
    };
    pub const Z: Direction = Direction {
        x: 0.,
        y: 0.,
        z: 1.,
    };

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

    pub fn to_cartesian(&self, length: Distance<f64>) -> CartesianCoordinates {
        CartesianCoordinates::new(self.x * length, self.y * length, self.z * length)
    }

    pub fn to_spherical(&self) -> SphericalCoordinates {
        SphericalCoordinates::cartesian_to_spherical((self.x, self.y, self.z))
    }

    pub fn x(&self) -> f64 {
        self.x
    }

    pub fn y(&self) -> f64 {
        self.y
    }

    pub fn z(&self) -> f64 {
        self.z
    }

    pub fn rotated(&self, angle: Angle<f64>, axis: &Direction) -> Direction {
        let (x, y, z) = rotated_tuple((self.x, self.y, self.z), angle, axis);
        Direction { x, y, z }
    }

    pub fn eq_within(&self, other: &Direction, accuracy: f64) -> bool {
        (self.x - other.x).abs() < accuracy
            && (self.y - other.y).abs() < accuracy
            && (self.z - other.z).abs() < accuracy
    }

    pub fn angle_to(&self, other: &Direction) -> Angle<f64> {
        let (ax, ay, az) = (self.x(), self.y(), self.z());
        let (bx, by, bz) = (other.x(), other.y(), other.z());

        let cosine_argument = ax * bx + ay * by + az * bz; //Directions have unit Distance
        if cosine_argument > 1. {
            //Saving acos from being called with an argument > 1 due to numerical instability
            return ANGLE_ZERO;
        } else if cosine_argument < -1. {
            return HALF_CIRC;
        }
        Angle::from_radians(cosine_argument.acos())
    }

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

    pub fn cross_product(&self, other: &Direction) -> Result<Direction, AstroCoordsError> {
        let (ax, ay, az) = (self.x, self.y, self.z);
        let (bx, by, bz) = (other.x(), other.y(), other.z());

        let cx = ay * bz - az * by;
        let cy = az * bx - ax * bz;
        let cz = ax * by - ay * bx;

        Direction::new(cx, cy, cz)
    }

    pub fn dot_product(&self, other: &Direction) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /*
     * This method is for example used to convert from equatorial coordinates to ecliptic coordinates.
     * It operates the following way:
     * 1. The vector is rotated around the X axis by the angle between new and old Z axis.
     * 2. The vector is rotated around old Z by the angle between new Z and old Y, projected on the old X-Y plane.
     * The result is the coordinates. The new X-axis still lies in the old X-Y plane.
     */
    pub fn active_rotation_to_new_z_axis(&self, new_z: &Direction) -> Direction {
        let angle_to_old_z = new_z.angle_to(&Self::Z);

        let axis_projected_onto_xy_plane = Direction::new(new_z.x(), new_z.y(), 0.);
        let mut polar_rotation_angle = ANGLE_ZERO;
        if let Ok(axis_projected_onto_xy_plane) = axis_projected_onto_xy_plane {
            polar_rotation_angle = axis_projected_onto_xy_plane.angle_to(&Self::Y);
            if axis_projected_onto_xy_plane.x() < 0. {
                polar_rotation_angle = -polar_rotation_angle;
            }
        }

        let mut dir = self.rotated(-angle_to_old_z, &Self::X);
        dir = dir.rotated(-polar_rotation_angle, &Self::Z);
        dir
    }

    pub fn passive_rotation_to_new_z_axis(&self, new_z: &Direction) -> Direction {
        let axis_projected_onto_xy_plane = Direction::new(new_z.x(), new_z.y(), 0.);
        let mut polar_rotation_angle = ANGLE_ZERO;
        if let Ok(axis_projected_onto_xy_plane) = axis_projected_onto_xy_plane {
            polar_rotation_angle = axis_projected_onto_xy_plane.angle_to(&Self::Y);
            if axis_projected_onto_xy_plane.x() < 0. {
                polar_rotation_angle = -polar_rotation_angle;
            }
        }

        let angle_to_old_z = new_z.angle_to(&Self::Z);

        let mut dir = self.rotated(polar_rotation_angle, &Self::Z);
        dir = dir.rotated(angle_to_old_z, &Self::X);
        dir
    }

    pub fn to_earth_equatorial(&self) -> EarthEquatorialCoordinates {
        let dir_in_equatorial = self.rotated(EARTH_AXIS_TILT, &Direction::X);
        let spherical = dir_in_equatorial.to_spherical();
        EarthEquatorialCoordinates::new(spherical.get_longitude(), spherical.get_latitude())
    }

    pub fn to_ecliptic(&self) -> EclipticCoordinates {
        EclipticCoordinates::new(self.to_spherical())
    }

    pub fn to_equatorial(&self, axis: Direction) -> EquatorialCoordinates {
        let dir_in_equatorial = self.passive_rotation_to_new_z_axis(&axis);
        let spherical = dir_in_equatorial.to_spherical();
        EquatorialCoordinates::new(spherical, axis)
    }
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
                    let cartesian = CartesianCoordinates::new(x, y, z);
                    let length = cartesian.length();
                    let expected_x = x / length;
                    let expected_y = y / length;
                    let expected_z = z / length;
                    let direction = cartesian.to_spherical().to_direction();

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
        let cartesian = CartesianCoordinates::new(x, y, z);
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
