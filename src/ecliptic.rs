use serde::{Deserialize, Serialize};
use simple_si_units::geometry::Angle;
use std::{fmt::Display, ops::Neg};

use super::{direction::Direction, spherical::SphericalCoordinates};

/*  The "absolute" reference we use for polar coordiantes is heliocentric ecliptic coordinates:
 * Longitude denotes the angle between the vernal equinox and the body, measured in the ecliptic plane.
 * Latitude denotes the angle between the ecliptic plane and the body.
 * https://en.wikipedia.org/wiki/Ecliptic_coordinate_system
 */
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EclipticCoordinates {
    pub(super) spherical: SphericalCoordinates,
}

impl EclipticCoordinates {
    pub const X_DIRECTION: EclipticCoordinates = EclipticCoordinates {
        spherical: SphericalCoordinates::X_DIRECTION,
    };

    pub const Y_DIRECTION: EclipticCoordinates = EclipticCoordinates {
        spherical: SphericalCoordinates::Y_DIRECTION,
    };

    pub const Z_DIRECTION: EclipticCoordinates = EclipticCoordinates {
        spherical: SphericalCoordinates::Z_DIRECTION,
    };

    pub const fn new(spherical: SphericalCoordinates) -> EclipticCoordinates {
        EclipticCoordinates { spherical }
    }

    pub fn normalize(&mut self) {
        self.spherical.normalize();
    }

    pub fn eq_within(&self, other: &EclipticCoordinates, accuracy: Angle<f64>) -> bool {
        self.spherical.eq_within(&other.spherical, accuracy)
    }

    pub fn get_longitude(&self) -> Angle<f64> {
        self.spherical.get_longitude()
    }

    pub fn get_latitude(&self) -> Angle<f64> {
        self.spherical.get_latitude()
    }

    pub fn get_spherical(&self) -> &SphericalCoordinates {
        &self.spherical
    }

    pub fn set_longitude(&mut self, longitude: Angle<f64>) {
        self.spherical.set_longitude(longitude);
    }

    pub fn set_latitude(&mut self, latitude: Angle<f64>) {
        self.spherical.set_latitude(latitude);
    }

    pub fn to_direction(&self) -> Direction {
        self.spherical.to_direction()
    }

    pub fn angle_to(&self, other: &Self) -> Angle<f64> {
        //TODO: There must be a more performant and stable way.
        self.to_direction().angle_to(&other.to_direction())
    }
}

impl Neg for &EclipticCoordinates {
    type Output = EclipticCoordinates;

    fn neg(self) -> EclipticCoordinates {
        EclipticCoordinates {
            spherical: -self.spherical,
        }
    }
}

impl Neg for EclipticCoordinates {
    type Output = EclipticCoordinates;

    fn neg(self) -> EclipticCoordinates {
        -&self
    }
}

impl Display for EclipticCoordinates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.spherical)
    }
}

#[cfg(test)]
mod tests {
    use simple_si_units::geometry::Angle;

    use crate::{angle_helper::angle_eq_within, direction::Direction};

    #[test]
    fn test_angle_function() {
        let numbers = [1., 0., -1., 3.];
        for x in numbers.iter() {
            for y in numbers.iter() {
                for z in numbers.iter() {
                    for angle in numbers.iter() {
                        let dir1 = Direction::new(*x, *y, *z);
                        if dir1.is_err() {
                            continue;
                        }
                        let angle = Angle::from_radians((*angle).abs());
                        let dir1 = dir1.unwrap();
                        let axis = dir1.some_orthogonal_vector();
                        let dir2 = dir1.rotated(angle, &axis);

                        let ecliptic1 = dir1.to_ecliptic();
                        let ecliptic2 = dir2.to_ecliptic();
                        let actual_angle = ecliptic1.angle_to(&ecliptic2);

                        println!("Expected: {}", angle);
                        println!("Actual: {}", actual_angle);
                        assert!(angle_eq_within(actual_angle, angle, Angle { rad: 1e-5 }));
                    }
                }
            }
        }
    }
}
