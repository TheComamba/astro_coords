//! This module contains the Spherical struct and its implementation.

use serde::{Deserialize, Serialize};
use simple_si_units::{base::Distance, geometry::Angle};
use std::{fmt::Display, ops::Neg};

use crate::{
    angle_helper::*, cartesian::Cartesian, earth_equatorial::EarthEquatorial,
    equatorial::Equatorial, error::AstroCoordsError, traits::*, NORMALIZATION_THRESHOLD,
};

use super::{
    direction::Direction,
    ecliptic::Ecliptic,
    ra_and_dec::{Declination, RightAscension, Sgn},
};

/// A struct representing coordinates on the surface of a sphere.
///
/// Contrary to the mathematical definition of spherical coordinates, the longitude is measured from the x-axis, and the latitude is measured from the xy-plane.
#[derive(Debug, Copy, Clone, PartialEq, Serialize, Deserialize)]
pub struct Spherical {
    /// The longitude of the Spherical struct measures the angle between the x-axis and the vector projected on the x-y plane.
    pub longitude: Angle<f64>,
    /// The latitude of the Spherical struct measures the angle between the vector and the x-y plane.
    pub latitude: Angle<f64>,
}

impl Spherical {
    /// The x-axis direction represented in spherical coordinates.
    pub const X_DIRECTION: Spherical = Spherical {
        longitude: ANGLE_ZERO,
        latitude: ANGLE_ZERO,
    };

    /// The y-axis direction represented in spherical coordinates.
    pub const Y_DIRECTION: Spherical = Spherical {
        longitude: QUARTER_CIRC,
        latitude: ANGLE_ZERO,
    };

    /// The z-axis direction represented in spherical coordinates.
    pub const Z_DIRECTION: Spherical = Spherical {
        longitude: ANGLE_ZERO,
        latitude: QUARTER_CIRC,
    };

    /// Creates a new Spherical struct from longitude and latitude input.
    ///
    /// Contrary to the mathematical definition of spherical coordinates, the longitude is measured from the x-axis, and the latitude is measured from the xy-plane.
    pub const fn new(longitude: Angle<f64>, latitude: Angle<f64>) -> Self {
        Self {
            longitude,
            latitude,
        }
    }

    /// Normalizes the longitude and latitude of the Spherical struct.
    ///
    /// The longitude is normalized to the range [-π, π), and the latitude is normalized to the range [-π/2, π/2].
    ///
    /// # Examples
    /// ```
    /// use astro_coords::spherical::Spherical;
    /// use simple_si_units::geometry::Angle;
    ///
    /// let mut coords = Spherical::new(Angle::from_degrees(190.), Angle::from_degrees(0.));
    /// coords.normalize();
    /// assert!((coords.longitude.to_degrees() + 170.).abs() < 1e-5);
    ///
    /// let mut coords = Spherical::new(Angle::from_degrees(10.), Angle::from_degrees(100.));
    /// coords.normalize();
    /// assert!((coords.longitude.to_degrees() + 170.).abs() < 1e-5);
    /// assert!((coords.latitude.to_degrees() - 80.).abs() < 1e-5);
    /// ```
    pub fn normalize(&mut self) {
        self.longitude = normalized_angle(self.longitude);
        self.latitude = normalized_angle(self.latitude);
        if self.latitude > QUARTER_CIRC {
            self.longitude += HALF_CIRC;
            self.longitude = normalized_angle(self.longitude);
            self.latitude = HALF_CIRC - self.latitude;
        } else if self.latitude < -QUARTER_CIRC {
            self.longitude += HALF_CIRC;
            self.longitude = normalized_angle(self.longitude);
            self.latitude = -HALF_CIRC - self.latitude;
        }
    }

    /// Rotates the Spherical struct around the specified axis by the specified angle.
    ///
    /// # Examples
    /// ```
    /// use astro_coords::spherical::Spherical;
    /// use astro_coords::direction::Direction;
    /// use simple_si_units::geometry::Angle;
    ///
    /// let mut coords = Spherical::X_DIRECTION;
    /// let rotated = coords.rotated(Angle::from_degrees(90.), &Direction::Z);
    /// assert!(rotated.eq_within(&Spherical::Y_DIRECTION, Angle::from_degrees(1e-5)));
    /// ```
    pub fn rotated(&self, angle: Angle<f64>, axis: &Direction) -> Spherical {
        self.to_direction().rotated(angle, axis).to_spherical()
    }

    /// Rotates the Spherical struct around the x-axis by the specified angle.
    ///
    /// This is an ever so tiny bit faster than `rotated()`.
    ///
    /// # Examples
    /// ```
    /// use astro_coords::spherical::Spherical;
    /// use simple_si_units::geometry::Angle;
    ///
    /// let mut coords = Spherical::Y_DIRECTION;
    /// let rotated = coords.rotated_x(Angle::from_degrees(90.));
    /// assert!((rotated.eq_within(&Spherical::Z_DIRECTION, Angle::from_degrees(1e-5))));
    /// ```
    pub fn rotated_x(&self, angle: Angle<f64>) -> Spherical {
        self.to_direction().rotated_x(angle).to_spherical()
    }

    /// Rotates the Spherical struct around the y-axis by the specified angle.
    ///
    /// This is an ever so tiny bit faster than `rotated()`.
    ///
    /// # Examples
    /// ```
    /// use astro_coords::spherical::Spherical;
    /// use simple_si_units::geometry::Angle;
    ///
    /// let mut coords = Spherical::Z_DIRECTION;
    /// let rotated = coords.rotated_y(Angle::from_degrees(90.));
    /// assert!((rotated.eq_within(&Spherical::X_DIRECTION, Angle::from_degrees(1e-5))));
    /// ```
    pub fn rotated_y(&self, angle: Angle<f64>) -> Spherical {
        self.to_direction().rotated_y(angle).to_spherical()
    }

    /// Rotates the Spherical struct around the z-axis by the specified angle.
    ///
    /// This is much faster than `rotated()`.
    ///
    /// # Examples
    /// ```
    /// use astro_coords::spherical::Spherical;
    /// use simple_si_units::geometry::Angle;
    ///
    /// let mut coords = Spherical::X_DIRECTION;
    /// let rotated = coords.rotated_z(Angle::from_degrees(90.));
    /// assert!((rotated.eq_within(&Spherical::Y_DIRECTION, Angle::from_degrees(1e-5))));
    /// ```
    pub fn rotated_z(&self, angle: Angle<f64>) -> Spherical {
        let mut rotated = self.clone();
        rotated.longitude += angle;
        rotated
    }

    /// Returns the spherical coordinates that result from actively rotating the spherical vector to the new z-axis, in a manner that preserves the old z-projection of the x-axis.
    ///
    /// This method is for example used to convert from equatorial coordinates to ecliptic coordinates.
    /// It operates in the following way:
    /// 1. The vector is rotated around the old x-axis by the angle between new and old z-axis.
    /// 2. The vector is rotated around the old z-axis by the angle between the new z-axis and the old y-axis, projected onto the old x-y plane.
    ///
    /// This is the inverse operation of `passive_rotation_to_new_z_axis`. See there for a somewhat intuitive example.
    pub fn active_rotation_to_new_z_axis(&self, new_z: &Self) -> Self {
        let (angle_to_old_z, polar_rotation_angle) =
            get_angle_to_old_z_and_polar_rotation_angle(new_z);
        self.rotated_x(-angle_to_old_z)
            .rotated_z(-polar_rotation_angle)
    }

    /// Returns the spherical coordinates that result from passively rotating the spherical vector to the new z-axis, in a manner that preserves the old z-projection of the x-axis.
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
    /// use astro_coords::spherical::Spherical;
    /// use simple_si_units::geometry::Angle;
    ///
    /// // Suppose it is summer solstice and the sun is in y-direction in the ecliptic coordinate system.
    /// let dir_of_sun_in_ecliptic = Spherical::Y_DIRECTION;
    ///
    /// // Now we want to express the sun's direction in earth equatorial coordinates.
    /// // The rotation axis of the earth expressed in ecliptic coordinates is given by:
    /// let earth_axis_tilt = Angle::from_degrees(23.44);
    /// let earth_rotation_axis_in_ecliptic = Spherical::Z_DIRECTION.rotated_x(-earth_axis_tilt);
    ///
    /// // The sun's direction in earth equatorial coordinates is then:
    /// let dir_of_sun_in_equatorial = dir_of_sun_in_ecliptic.passive_rotation_to_new_z_axis(&earth_rotation_axis_in_ecliptic);
    ///
    /// // At summer solstice, the sun is highest in the sky in the northern hemisphere, so its longitude is a quarter turn, and its latitude is positive.
    /// println!("{}", dir_of_sun_in_equatorial);
    /// assert!((dir_of_sun_in_equatorial.longitude.to_degrees() - 90.).abs() < 1e-5);
    /// assert!(dir_of_sun_in_equatorial.latitude.to_degrees() > 0.);
    /// ```
    pub fn passive_rotation_to_new_z_axis(&self, new_z: &Self) -> Self {
        let (angle_to_old_z, polar_rotation_angle) =
            get_angle_to_old_z_and_polar_rotation_angle(new_z);
        self.rotated_z(polar_rotation_angle)
            .rotated_x(angle_to_old_z)
    }

    /// Checks if the longitude and latitude of the Spherical struct are quivalent to those of another Spherical struct within a certain accuracy.
    ///
    /// # Examples
    /// ```
    /// use astro_coords::spherical::Spherical;
    /// use simple_si_units::geometry::Angle;
    ///
    /// let mut c1 = Spherical::new(Angle::from_degrees(10.), Angle::from_degrees(100.));
    /// let mut c2 = Spherical::new(Angle::from_degrees(-170.), Angle::from_degrees(80.));
    /// assert!(c1.eq_within(&c2, Angle::from_degrees(1.)));
    /// ```
    pub fn eq_within(&self, other: &Self, accuracy: Angle<f64>) -> bool {
        let northpole_latitude = QUARTER_CIRC;
        let southpole_latitude = -QUARTER_CIRC;
        let mut copy = *self;
        let mut other_copy = *other;
        copy.normalize();
        other_copy.normalize();
        let latitudes_equal = angle_eq_within(copy.latitude, other_copy.latitude, accuracy);
        let is_pole = angle_eq_within(copy.latitude, northpole_latitude, accuracy)
            || angle_eq_within(copy.latitude, southpole_latitude, accuracy);
        let longitudes_equal = if is_pole {
            true
        } else {
            angle_eq_within(copy.longitude, other_copy.longitude, accuracy)
        };
        latitudes_equal && longitudes_equal
    }

    /// Calculates the angle between two points on a sphere (great-circle distance) given in spherical coordinates.
    ///
    /// This method uses the Haversine formula:
    /// https://en.wikipedia.org/wiki/Haversine_formula
    ///
    /// # Example
    /// ```
    /// use astro_coords::spherical::Spherical;
    /// use simple_si_units::geometry::Angle;
    ///
    /// let a = Spherical::new(Angle::from_degrees(0.), Angle::from_degrees(0.));
    /// let b = Spherical::new(Angle::from_degrees(90.), Angle::from_degrees(0.));
    /// let angle = a.angle_to(&b);
    /// assert!((angle - Angle::from_degrees(90.)).rad.abs() < 1e-5);
    /// ```
    pub fn angle_to(&self, other: &Self) -> Angle<f64> {
        let lat_diff = other.latitude - self.latitude;
        let lat_sum = other.latitude + self.latitude;
        let long_diff = other.longitude - self.longitude;
        let hav_lat_diff = haversine(lat_diff);
        let hav_lat_sum = haversine(lat_sum);
        let hav_long_diff = haversine(long_diff);
        let arg = hav_lat_diff + (1. - hav_lat_diff - hav_lat_sum) * hav_long_diff;
        arcushaversine(arg)
    }

    pub(crate) fn cartesian_to_spherical(cart: (f64, f64, f64)) -> Result<Self, AstroCoordsError> {
        let (x, y, z) = cart;
        if x * x + y * y + z * z < NORMALIZATION_THRESHOLD {
            return Err(AstroCoordsError::NormalizingZeroVector);
        }

        let longitude = Angle::from_radians(y.atan2(x));
        let latitude = Angle::from_radians(z.atan2((x * x + y * y).sqrt()));
        Ok(Self {
            longitude,
            latitude,
        })
    }

    /// Constructs Cartesian with the provided length pointing in the direction specified by the Spherical struct.
    pub fn to_cartesian(&self, length: Distance<f64>) -> Cartesian {
        self.to_direction().to_cartesian(length)
    }

    /// Constructs Direction with the same direction as the Spherical struct.
    ///
    /// # Examples
    /// ```
    /// use astro_coords::spherical::Spherical;
    /// use astro_coords::direction::Direction;
    ///
    /// let x = Spherical::X_DIRECTION.to_direction();
    /// let y = Spherical::Y_DIRECTION.to_direction();
    /// let z = Spherical::Z_DIRECTION.to_direction();
    /// assert!(x.eq_within(&Direction::X, 1e-5));
    /// assert!(y.eq_within(&Direction::Y, 1e-5));
    /// assert!(z.eq_within(&Direction::Z, 1e-5));
    /// ```
    pub fn to_direction(&self) -> Direction {
        let x = self.longitude.rad.cos() * self.latitude.rad.cos();
        let y = self.longitude.rad.sin() * self.latitude.rad.cos();
        let z = self.latitude.rad.sin();
        Direction { x, y, z }
    }

    pub fn to_earth_equatorial(&self) -> EarthEquatorial {
        self.to_direction().to_earth_equatorial()
    }

    pub fn to_ecliptic(&self) -> Ecliptic {
        Ecliptic { spherical: *self }
    }

    pub fn to_equatorial(&self, axis: Direction) -> Equatorial {
        self.to_direction().to_equatorial(axis)
    }

    pub fn to_ra_and_dec(&self) -> (RightAscension, Declination) {
        let mut ra_remainder = normalized_angle(self.longitude).to_degrees();
        if ra_remainder < 0. {
            ra_remainder += 360.;
        }
        let ra_hours = (ra_remainder / 15.).floor() as u8;
        ra_remainder -= ra_hours as f64 * 15.;
        let ra_minutes = (ra_remainder / 15. * 60.).floor() as u8;
        ra_remainder -= ra_minutes as f64 / 60. * 15.;
        let ra_seconds = (ra_remainder / 15. * 3600.).floor();
        let ra = RightAscension::new(ra_hours, ra_minutes, ra_seconds);

        let mut dec_remainder = self.latitude.to_degrees();
        let sign = if dec_remainder < 0. {
            dec_remainder = dec_remainder.abs();
            Sgn::Neg
        } else {
            Sgn::Pos
        };
        let dec_degrees = dec_remainder.floor() as u8;
        dec_remainder -= dec_degrees as f64;
        let dec_minutes = (dec_remainder * 60.).floor() as u8;
        dec_remainder -= dec_minutes as f64 / 60.;
        let dec_seconds = (dec_remainder * 3600.).floor();
        let dec = Declination::new(sign, dec_degrees, dec_minutes, dec_seconds);

        (ra, dec)
    }
}

fn get_angle_to_old_z_and_polar_rotation_angle(new_z: &Spherical) -> (Angle<f64>, Angle<f64>) {
    let angle_to_old_z = QUARTER_CIRC - new_z.latitude;
    let is_polar =
        angle_to_old_z.rad.abs() < 1e-20 || (angle_to_old_z - HALF_CIRC).rad.abs() < 1e-20;
    let polar_rotation_angle = if is_polar {
        Angle { rad: 0. }
    } else {
        QUARTER_CIRC - new_z.longitude
    };
    (angle_to_old_z, polar_rotation_angle)
}

impl Neg for &Spherical {
    type Output = Spherical;

    fn neg(self) -> Spherical {
        let mut longitude = self.longitude + HALF_CIRC;
        longitude = normalized_angle(longitude);
        let latitude = -self.latitude;
        Spherical {
            longitude,
            latitude,
        }
    }
}

impl Neg for Spherical {
    type Output = Spherical;

    fn neg(self) -> Spherical {
        -&self
    }
}

impl Display for Spherical {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({} long, {} lat)", self.longitude, self.latitude)
    }
}

#[cfg(test)]
mod tests {
    use simple_si_units::base::Distance;
    use std::f64::consts::PI;

    use crate::cartesian::Cartesian;

    use super::*;

    const TEST_ACCURACY: f64 = 1e-5;
    const ANGLE_TEST_ACCURACY: Angle<f64> = Angle { rad: TEST_ACCURACY };

    #[test]
    fn test_from_cartesian() {
        let cartesian = Cartesian::new(
            Distance::from_meters(1.),
            Distance::from_meters(1.),
            Distance::from_meters(1.),
        );
        let expected: Spherical = Spherical {
            longitude: Angle::from_degrees(45.),
            latitude: Angle::from_degrees(90. - 54.7356),
        };
        let actual = cartesian.to_spherical().unwrap();
        println!("{}, expected: {}, actual: {}", cartesian, expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = Cartesian::new(
            Distance::from_meters(1.),
            Distance::from_meters(1.),
            Distance::from_meters(-1.),
        );
        let expected = Spherical {
            longitude: Angle::from_degrees(45.),
            latitude: Angle::from_degrees(-90. + 54.7356),
        };
        let actual = cartesian.to_spherical().unwrap();
        println!("{}, expected: {}, actual: {}", cartesian, expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = Cartesian::new(
            Distance::from_meters(1.),
            Distance::from_meters(-1.),
            Distance::from_meters(1.),
        );
        let expected = Spherical {
            longitude: Angle::from_degrees(-45.),
            latitude: Angle::from_degrees(90. - 54.7356),
        };
        let actual = cartesian.to_spherical().unwrap();
        println!("{}, expected: {}, actual: {}", cartesian, expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = Cartesian::new(
            Distance::from_meters(1.),
            Distance::from_meters(-1.),
            Distance::from_meters(-1.),
        );
        let expected = Spherical {
            longitude: Angle::from_degrees(-45.),
            latitude: Angle::from_degrees(-90. + 54.7356),
        };
        let actual = cartesian.to_spherical().unwrap();
        println!("{}, expected: {}, actual: {}", cartesian, expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = Cartesian::new(
            Distance::from_meters(-1.),
            Distance::from_meters(1.),
            Distance::from_meters(1.),
        );
        let expected = Spherical {
            longitude: Angle::from_degrees(135.),
            latitude: Angle::from_degrees(90. - 54.7356),
        };
        let actual = cartesian.to_spherical().unwrap();
        println!("{}, expected: {}, actual: {}", cartesian, expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = Cartesian::new(
            Distance::from_meters(-1.),
            Distance::from_meters(1.),
            Distance::from_meters(-1.),
        );
        let expected = Spherical {
            longitude: Angle::from_degrees(135.),
            latitude: Angle::from_degrees(-90. + 54.7356),
        };
        let actual = cartesian.to_spherical().unwrap();
        println!("{}, expected: {}, actual: {}", cartesian, expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = Cartesian::new(
            Distance::from_meters(-1.),
            Distance::from_meters(-1.),
            Distance::from_meters(1.),
        );
        let expected = Spherical {
            longitude: Angle::from_degrees(-135.),
            latitude: Angle::from_degrees(90. - 54.7356),
        };
        let actual = cartesian.to_spherical().unwrap();
        println!("{}, expected: {}, actual: {}", cartesian, expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = Cartesian::new(
            Distance::from_meters(-1.),
            Distance::from_meters(-1.),
            Distance::from_meters(-1.),
        );
        let expected = Spherical {
            longitude: Angle::from_degrees(-135.),
            latitude: Angle::from_degrees(-90. + 54.7356),
        };
        let actual = cartesian.to_spherical().unwrap();
        println!("{}, expected: {}, actual: {}", cartesian, expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    #[test]
    fn test_unary_minus() {
        let x = Spherical::X_DIRECTION;
        let y = Spherical::Y_DIRECTION;
        let z = Spherical::Z_DIRECTION;
        let xyz = Spherical {
            longitude: Angle::from_radians(PI / 4.),
            latitude: Angle::from_radians(PI / 4.),
        };
        let expected_minus_x = Spherical {
            longitude: HALF_CIRC,
            latitude: ANGLE_ZERO,
        };
        let expected_minus_y = Spherical {
            longitude: -QUARTER_CIRC,
            latitude: ANGLE_ZERO,
        };
        let expected_minus_z = Spherical {
            longitude: ANGLE_ZERO,
            latitude: -QUARTER_CIRC,
        };
        let expected_minus_xyz = Spherical {
            longitude: Angle::from_radians(-PI * 3. / 4.),
            latitude: Angle::from_radians(-PI / 4.),
        };

        let minus_x = -x;
        println!("-x, expected: {}, actual: {}", expected_minus_x, minus_x);
        assert!(minus_x.eq_within(&expected_minus_x, ANGLE_TEST_ACCURACY));

        let minus_y = -y;
        println!("-y, expected: {}, actual: {}", expected_minus_y, minus_y);
        assert!(minus_y.eq_within(&expected_minus_y, ANGLE_TEST_ACCURACY));

        let minus_z = -z;
        println!("-z, expected: {}, actual: {}", expected_minus_z, minus_z);
        assert!(minus_z.eq_within(&expected_minus_z, ANGLE_TEST_ACCURACY));

        let minus_xyz = -xyz;
        println!(
            "-xyz, expected: {}, actual: {}",
            expected_minus_xyz, minus_xyz
        );
        assert!(minus_xyz.eq_within(&expected_minus_xyz, ANGLE_TEST_ACCURACY));
    }

    #[test]
    fn test_eq_within() {
        const LOCAL_TEST_ANGLE_ACCURACY: Angle<f64> = Angle {
            rad: 10. * TEST_ACCURACY,
        };

        let small_offsets = vec![
            -LOCAL_TEST_ANGLE_ACCURACY / 100.,
            ANGLE_ZERO,
            LOCAL_TEST_ANGLE_ACCURACY / 100.,
        ];
        let large_offsets = vec![-FULL_CIRC, ANGLE_ZERO, FULL_CIRC, 100. * FULL_CIRC];
        let directions = vec![
            Spherical::X_DIRECTION,
            Spherical::Y_DIRECTION,
            Spherical::Z_DIRECTION,
            -Spherical::X_DIRECTION,
            -Spherical::Y_DIRECTION,
            -Spherical::Z_DIRECTION,
        ];
        for direction in directions.clone() {
            for small_offset in small_offsets.clone() {
                for large_offset in large_offsets.clone() {
                    let longitude1 = direction.longitude;
                    let latitude1 = direction.latitude;
                    let longitude2 = longitude1 + large_offset + small_offset;
                    let latitude2 = latitude1 + large_offset + small_offset;
                    let coords1 = Spherical {
                        longitude: longitude2,
                        latitude: latitude1,
                    };
                    let coords2 = Spherical {
                        longitude: longitude1,
                        latitude: latitude2,
                    };
                    let coords3 = Spherical {
                        longitude: longitude2,
                        latitude: latitude2,
                    };
                    println!("Expecting {} == {}", direction, coords1);
                    assert![direction.eq_within(&coords1, LOCAL_TEST_ANGLE_ACCURACY)];
                    println!("Expecting {} == {}", direction, coords2);
                    assert![direction.eq_within(&coords2, LOCAL_TEST_ANGLE_ACCURACY)];
                    println!("Expecting {} == {}", direction, coords3);
                    assert![direction.eq_within(&coords3, LOCAL_TEST_ANGLE_ACCURACY)];
                }
            }
        }
    }

    #[test]
    fn test_normalization_two_pi_offsets() {
        const LOCAL_TEST_ANGLE_ACCURACY: Angle<f64> = Angle {
            rad: 10. * TEST_ACCURACY,
        };

        let longitudes = vec![
            -0.75 * PI,
            -0.5 * PI,
            -0.25 * PI,
            0.,
            0.25 * PI,
            0.5 * PI,
            0.75 * PI,
            PI,
            1.25 * PI,
            1.5 * PI,
            1.75 * PI,
            2. * PI,
        ];
        let latitudes = vec![-0.25 * PI, 0., 0.25 * PI];
        let offsets = vec![-FULL_CIRC, ANGLE_ZERO, FULL_CIRC, 100. * FULL_CIRC];
        for longitude in longitudes.clone() {
            for latitude in latitudes.clone() {
                for offset in offsets.clone() {
                    let longitude1 = Angle::from_radians(longitude);
                    let latitude1 = Angle::from_radians(latitude);
                    let longitude2 = longitude1 + offset;
                    let latitude2 = latitude1 + offset;
                    let mut coords1 = Spherical {
                        longitude: longitude1,
                        latitude: latitude1,
                    };
                    let mut coords2 = Spherical {
                        longitude: longitude1,
                        latitude: latitude2,
                    };
                    let mut coords3 = Spherical {
                        longitude: longitude2,
                        latitude: latitude1,
                    };
                    let mut coords4 = Spherical {
                        longitude: longitude2,
                        latitude: latitude2,
                    };
                    coords1.normalize();
                    coords2.normalize();
                    coords3.normalize();
                    coords4.normalize();
                    assert!(coords1.eq_within(&coords2, LOCAL_TEST_ANGLE_ACCURACY));
                    assert!(coords1.eq_within(&coords3, LOCAL_TEST_ANGLE_ACCURACY));
                    assert!(coords1.eq_within(&coords4, LOCAL_TEST_ANGLE_ACCURACY));
                }
            }
        }
    }

    #[test]
    fn test_normalization_crossing_poles() {
        let mut coord = Spherical {
            longitude: ANGLE_ZERO,
            latitude: Angle::from_radians(3. / 4. * PI),
        };
        let expected = Spherical {
            longitude: HALF_CIRC,
            latitude: Angle::from_radians(PI / 4.),
        };
        coord.normalize();
        println!("expected: {}, actual: {}", expected, coord);
        assert!(coord.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let mut coord = Spherical {
            longitude: ANGLE_ZERO,
            latitude: Angle::from_radians(-3. / 4. * PI),
        };
        let expected = Spherical {
            longitude: HALF_CIRC,
            latitude: Angle::from_radians(-PI / 4.),
        };
        coord.normalize();
        println!("expected: {}, actual: {}", expected, coord);
        assert!(coord.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let mut coord = Spherical {
            longitude: HALF_CIRC,
            latitude: Angle::from_radians(3. / 4. * PI),
        };
        let expected = Spherical {
            longitude: ANGLE_ZERO,
            latitude: Angle::from_radians(PI / 4.),
        };
        coord.normalize();
        println!("expected: {}, actual: {}", expected, coord);
        assert!(coord.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let mut coord = Spherical {
            longitude: HALF_CIRC,
            latitude: Angle::from_radians(-3. / 4. * PI),
        };
        let expected = Spherical {
            longitude: ANGLE_ZERO,
            latitude: Angle::from_radians(-PI / 4.),
        };
        coord.normalize();
        println!("expected: {}, actual: {}", expected, coord);
        assert!(coord.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let mut coord = Spherical {
            longitude: QUARTER_CIRC,
            latitude: Angle::from_radians(3. / 4. * PI),
        };
        let expected = Spherical {
            longitude: Angle::from_radians(3. / 2. * PI),
            latitude: Angle::from_radians(PI / 4.),
        };
        coord.normalize();
        println!("expected: {}, actual: {}", expected, coord);
        assert!(coord.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let mut coord = Spherical {
            longitude: QUARTER_CIRC,
            latitude: Angle::from_radians(-3. / 4. * PI),
        };
        let expected = Spherical {
            longitude: Angle::from_radians(3. / 2. * PI),
            latitude: Angle::from_radians(-PI / 4.),
        };
        coord.normalize();
        println!("expected: {}, actual: {}", expected, coord);
        assert!(coord.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let mut coord = Spherical {
            longitude: Angle::from_radians(3. / 2. * PI),
            latitude: Angle::from_radians(3. / 4. * PI),
        };
        let expected = Spherical {
            longitude: QUARTER_CIRC,
            latitude: Angle::from_radians(PI / 4.),
        };
        coord.normalize();
        println!("expected: {}, actual: {}", expected, coord);
        assert!(coord.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    #[test]
    fn roundtrips_to_ra_and_dec() {
        let local_test_accuracy = 10. * ANGLE_TEST_ACCURACY;
        const STEP: usize = 100;
        for i in 0..STEP {
            for j in 0..STEP {
                let ra_angle = Angle::from_radians(2. * PI * i as f64 / STEP as f64);
                let dec_angle = Angle::from_radians(PI * j as f64 / STEP as f64 - PI / 2.);
                let spherical = Spherical {
                    longitude: ra_angle,
                    latitude: dec_angle,
                };
                let (ra, dec) = spherical.to_ra_and_dec();
                let spherical2 = Spherical::new(ra.to_angle(), dec.to_angle());
                println!("spherical: {}, spherical2: {}", spherical, spherical2);
                assert!(spherical.eq_within(&spherical2, local_test_accuracy));
            }
        }
    }
}
