use serde::{Deserialize, Serialize};
use simple_si_units::geometry::Angle;
use std::ops::Neg;

use crate::angle_helper::*;

use super::{
    declination::{Declination, Sgn},
    direction::Direction,
    ecliptic::EclipticCoordinates,
    right_ascension::RightAscension,
};

#[derive(Debug, Copy, Clone, PartialEq, Serialize, Deserialize)]
pub struct SphericalCoordinates {
    longitude: Angle<f64>,
    latitude: Angle<f64>,
}

impl SphericalCoordinates {
    pub const X_DIRECTION: SphericalCoordinates = SphericalCoordinates {
        longitude: ANGLE_ZERO,
        latitude: ANGLE_ZERO,
    };

    pub const Y_DIRECTION: SphericalCoordinates = SphericalCoordinates {
        longitude: QUARTER_CIRC,
        latitude: ANGLE_ZERO,
    };

    pub const Z_DIRECTION: SphericalCoordinates = SphericalCoordinates {
        longitude: ANGLE_ZERO,
        latitude: QUARTER_CIRC,
    };

    pub const fn new(longitude: Angle<f64>, latitude: Angle<f64>) -> Self {
        Self {
            longitude,
            latitude,
        }
    }
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

    pub(crate) fn eq_within(&self, other: &Self, accuracy: Angle<f64>) -> bool {
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

    pub fn get_longitude(&self) -> Angle<f64> {
        self.longitude
    }

    pub fn get_latitude(&self) -> Angle<f64> {
        self.latitude
    }

    pub fn set_longitude(&mut self, longitude: Angle<f64>) {
        self.longitude = longitude;
    }

    pub fn set_latitude(&mut self, latitude: Angle<f64>) {
        self.latitude = latitude;
    }

    pub(super) fn cartesian_to_spherical(cart: (f64, f64, f64)) -> Self {
        let (x, y, z) = cart;
        let longitude = Angle::from_radians(y.atan2(x));
        let latitude = Angle::from_radians(z.atan2((x * x + y * y).sqrt()));
        Self {
            longitude,
            latitude,
        }
    }

    pub fn to_direction(&self) -> Direction {
        let x = self.get_longitude().rad.cos() * self.get_latitude().rad.cos();
        let y = self.get_longitude().rad.sin() * self.get_latitude().rad.cos();
        let z = self.get_latitude().rad.sin();
        Direction { x, y, z }
    }

    pub fn to_ecliptic(&self) -> EclipticCoordinates {
        EclipticCoordinates { spherical: *self }
    }

    pub fn to_ra_and_dec(&self) -> (RightAscension, Declination) {
        let mut ra_remainder = self.longitude.to_degrees();
        let ra_hours = (ra_remainder / 15.).floor() as i8;
        ra_remainder -= ra_hours as f64 * 15.;
        let ra_minutes = (ra_remainder / 15. * 60.).floor() as i8;
        ra_remainder -= ra_minutes as f64 / 60. * 15.;
        let ra_seconds = (ra_remainder / 15. * 3600.).floor() as i8;
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
        let dec_seconds = (dec_remainder * 3600.).floor() as u8;
        let dec = Declination::new(sign, dec_degrees, dec_minutes, dec_seconds);

        (ra, dec)
    }
}

impl Neg for &SphericalCoordinates {
    type Output = SphericalCoordinates;

    fn neg(self) -> SphericalCoordinates {
        let mut longitude = self.longitude + HALF_CIRC;
        longitude = normalized_angle(longitude);
        let latitude = -self.latitude;
        SphericalCoordinates {
            longitude,
            latitude,
        }
    }
}

impl Neg for SphericalCoordinates {
    type Output = SphericalCoordinates;

    fn neg(self) -> SphericalCoordinates {
        -&self
    }
}

#[cfg(test)]
mod tests {
    use simple_si_units::base::Distance;
    use std::f64::consts::PI;

    use crate::cartesian::CartesianCoordinates;

    use super::*;

    const TEST_ACCURACY: f64 = 1e-5;
    const ANGLE_TEST_ACCURACY: Angle<f64> = Angle { rad: TEST_ACCURACY };

    #[test]
    fn test_from_cartesian() {
        let cartesian = CartesianCoordinates::new(
            Distance::from_meters(1.),
            Distance::from_meters(1.),
            Distance::from_meters(1.),
        );
        let expected: SphericalCoordinates = SphericalCoordinates {
            longitude: Angle::from_degrees(45.),
            latitude: Angle::from_degrees(90. - 54.7356),
        };
        let actual = cartesian.to_spherical();
        println!(
            "{}, expected: {}, actual: {}",
            cartesian.astro_display(),
            expected,
            actual
        );
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = CartesianCoordinates::new(
            Distance::from_meters(1.),
            Distance::from_meters(1.),
            Distance::from_meters(-1.),
        );
        let expected = SphericalCoordinates {
            longitude: Angle::from_degrees(45.),
            latitude: Angle::from_degrees(-90. + 54.7356),
        };
        let actual = cartesian.to_spherical();
        println!(
            "{}, expected: {}, actual: {}",
            cartesian.astro_display(),
            expected,
            actual
        );
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = CartesianCoordinates::new(
            Distance::from_meters(1.),
            Distance::from_meters(-1.),
            Distance::from_meters(1.),
        );
        let expected = SphericalCoordinates {
            longitude: Angle::from_degrees(-45.),
            latitude: Angle::from_degrees(90. - 54.7356),
        };
        let actual = cartesian.to_spherical();
        println!(
            "{}, expected: {}, actual: {}",
            cartesian.astro_display(),
            expected,
            actual
        );
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = CartesianCoordinates::new(
            Distance::from_meters(1.),
            Distance::from_meters(-1.),
            Distance::from_meters(-1.),
        );
        let expected = SphericalCoordinates {
            longitude: Angle::from_degrees(-45.),
            latitude: Angle::from_degrees(-90. + 54.7356),
        };
        let actual = cartesian.to_spherical();
        println!(
            "{}, expected: {}, actual: {}",
            cartesian.astro_display(),
            expected,
            actual
        );
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = CartesianCoordinates::new(
            Distance::from_meters(-1.),
            Distance::from_meters(1.),
            Distance::from_meters(1.),
        );
        let expected = SphericalCoordinates {
            longitude: Angle::from_degrees(135.),
            latitude: Angle::from_degrees(90. - 54.7356),
        };
        let actual = cartesian.to_spherical();
        println!(
            "{}, expected: {}, actual: {}",
            cartesian.astro_display(),
            expected,
            actual
        );
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = CartesianCoordinates::new(
            Distance::from_meters(-1.),
            Distance::from_meters(1.),
            Distance::from_meters(-1.),
        );
        let expected = SphericalCoordinates {
            longitude: Angle::from_degrees(135.),
            latitude: Angle::from_degrees(-90. + 54.7356),
        };
        let actual = cartesian.to_spherical();
        println!(
            "{}, expected: {}, actual: {}",
            cartesian.astro_display(),
            expected,
            actual
        );
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = CartesianCoordinates::new(
            Distance::from_meters(-1.),
            Distance::from_meters(-1.),
            Distance::from_meters(1.),
        );
        let expected = SphericalCoordinates {
            longitude: Angle::from_degrees(-135.),
            latitude: Angle::from_degrees(90. - 54.7356),
        };
        let actual = cartesian.to_spherical();
        println!(
            "{}, expected: {}, actual: {}",
            cartesian.astro_display(),
            expected,
            actual
        );
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let cartesian = CartesianCoordinates::new(
            Distance::from_meters(-1.),
            Distance::from_meters(-1.),
            Distance::from_meters(-1.),
        );
        let expected = SphericalCoordinates {
            longitude: Angle::from_degrees(-135.),
            latitude: Angle::from_degrees(-90. + 54.7356),
        };
        let actual = cartesian.to_spherical();
        println!(
            "{}, expected: {}, actual: {}",
            cartesian.astro_display(),
            expected,
            actual
        );
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    #[test]
    fn test_unary_minus() {
        let x = SphericalCoordinates::X_DIRECTION;
        let y = SphericalCoordinates::Y_DIRECTION;
        let z = SphericalCoordinates::Z_DIRECTION;
        let xyz = SphericalCoordinates {
            longitude: Angle::from_radians(PI / 4.),
            latitude: Angle::from_radians(PI / 4.),
        };
        let expected_minus_x = SphericalCoordinates {
            longitude: HALF_CIRC,
            latitude: ANGLE_ZERO,
        };
        let expected_minus_y = SphericalCoordinates {
            longitude: -QUARTER_CIRC,
            latitude: ANGLE_ZERO,
        };
        let expected_minus_z = SphericalCoordinates {
            longitude: ANGLE_ZERO,
            latitude: -QUARTER_CIRC,
        };
        let expected_minus_xyz = SphericalCoordinates {
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
            SphericalCoordinates::X_DIRECTION,
            SphericalCoordinates::Y_DIRECTION,
            SphericalCoordinates::Z_DIRECTION,
            -SphericalCoordinates::X_DIRECTION,
            -SphericalCoordinates::Y_DIRECTION,
            -SphericalCoordinates::Z_DIRECTION,
        ];
        for direction in directions.clone() {
            for small_offset in small_offsets.clone() {
                for large_offset in large_offsets.clone() {
                    let longitude1 = direction.longitude;
                    let latitude1 = direction.latitude;
                    let longitude2 = longitude1 + large_offset + small_offset;
                    let latitude2 = latitude1 + large_offset + small_offset;
                    let coords1 = SphericalCoordinates {
                        longitude: longitude2,
                        latitude: latitude1,
                    };
                    let coords2 = SphericalCoordinates {
                        longitude: longitude1,
                        latitude: latitude2,
                    };
                    let coords3 = SphericalCoordinates {
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
                    let mut coords1 = SphericalCoordinates {
                        longitude: longitude1,
                        latitude: latitude1,
                    };
                    let mut coords2 = SphericalCoordinates {
                        longitude: longitude1,
                        latitude: latitude2,
                    };
                    let mut coords3 = SphericalCoordinates {
                        longitude: longitude2,
                        latitude: latitude1,
                    };
                    let mut coords4 = SphericalCoordinates {
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
        let mut coord = SphericalCoordinates {
            longitude: ANGLE_ZERO,
            latitude: Angle::from_radians(3. / 4. * PI),
        };
        let expected = SphericalCoordinates {
            longitude: HALF_CIRC,
            latitude: Angle::from_radians(PI / 4.),
        };
        coord.normalize();
        println!("expected: {}, actual: {}", expected, coord);
        assert!(coord.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let mut coord = SphericalCoordinates {
            longitude: ANGLE_ZERO,
            latitude: Angle::from_radians(-3. / 4. * PI),
        };
        let expected = SphericalCoordinates {
            longitude: HALF_CIRC,
            latitude: Angle::from_radians(-PI / 4.),
        };
        coord.normalize();
        println!("expected: {}, actual: {}", expected, coord);
        assert!(coord.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let mut coord = SphericalCoordinates {
            longitude: HALF_CIRC,
            latitude: Angle::from_radians(3. / 4. * PI),
        };
        let expected = SphericalCoordinates {
            longitude: ANGLE_ZERO,
            latitude: Angle::from_radians(PI / 4.),
        };
        coord.normalize();
        println!("expected: {}, actual: {}", expected, coord);
        assert!(coord.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let mut coord = SphericalCoordinates {
            longitude: HALF_CIRC,
            latitude: Angle::from_radians(-3. / 4. * PI),
        };
        let expected = SphericalCoordinates {
            longitude: ANGLE_ZERO,
            latitude: Angle::from_radians(-PI / 4.),
        };
        coord.normalize();
        println!("expected: {}, actual: {}", expected, coord);
        assert!(coord.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let mut coord = SphericalCoordinates {
            longitude: QUARTER_CIRC,
            latitude: Angle::from_radians(3. / 4. * PI),
        };
        let expected = SphericalCoordinates {
            longitude: Angle::from_radians(3. / 2. * PI),
            latitude: Angle::from_radians(PI / 4.),
        };
        coord.normalize();
        println!("expected: {}, actual: {}", expected, coord);
        assert!(coord.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let mut coord = SphericalCoordinates {
            longitude: QUARTER_CIRC,
            latitude: Angle::from_radians(-3. / 4. * PI),
        };
        let expected = SphericalCoordinates {
            longitude: Angle::from_radians(3. / 2. * PI),
            latitude: Angle::from_radians(-PI / 4.),
        };
        coord.normalize();
        println!("expected: {}, actual: {}", expected, coord);
        assert!(coord.eq_within(&expected, ANGLE_TEST_ACCURACY));

        let mut coord = SphericalCoordinates {
            longitude: Angle::from_radians(3. / 2. * PI),
            latitude: Angle::from_radians(3. / 4. * PI),
        };
        let expected = SphericalCoordinates {
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
                let spherical = SphericalCoordinates {
                    longitude: ra_angle,
                    latitude: dec_angle,
                };
                let (ra, dec) = spherical.to_ra_and_dec();
                let spherical2 = SphericalCoordinates::new(ra.to_angle(), dec.to_angle());
                println!("spherical: {}, spherical2: {}", spherical, spherical2);
                assert!(spherical.eq_within(&spherical2, local_test_accuracy));
            }
        }
    }
}
