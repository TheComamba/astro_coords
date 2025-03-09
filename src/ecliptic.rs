//! This module contains the Ecliptic struct and its implementation.

use std::{fmt::Display, ops::Neg};

use serde::{Deserialize, Serialize};
use uom::si::f64::{Angle, Length};

use crate::{cartesian::Cartesian, equatorial::Equatorial};

use super::{direction::Direction, spherical::Spherical};

/// Ecliptic coordinates are a spherical coordinate system that is centered on the sun.
///
/// The "absolute" reference we use for polar coordiantes is heliocentric ecliptic coordinates:
/// Longitude denotes the angle between the vernal equinox and the body, measured in the ecliptic plane.
/// Latitude denotes the angle between the ecliptic plane and the body.
/// Compare https://en.wikipedia.org/wiki/Ecliptic_coordinate_system
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Ecliptic {
    pub spherical: Spherical,
}

impl Ecliptic {
    pub const X_DIRECTION: Ecliptic = Ecliptic {
        spherical: Spherical::X_DIRECTION,
    };

    pub const Y_DIRECTION: Ecliptic = Ecliptic {
        spherical: Spherical::Y_DIRECTION,
    };

    pub const Z_DIRECTION: Ecliptic = Ecliptic {
        spherical: Spherical::Z_DIRECTION,
    };

    pub const fn new(spherical: Spherical) -> Ecliptic {
        Ecliptic { spherical }
    }

    pub fn normalize(&mut self) {
        self.spherical.normalize();
    }

    pub fn eq_within(&self, other: &Ecliptic, accuracy: Angle) -> bool {
        self.spherical.eq_within(&other.spherical, accuracy)
    }

    pub fn to_cartesian(&self, length: Length) -> Cartesian {
        self.to_direction().to_cartesian(length)
    }

    pub fn to_direction(&self) -> Direction {
        self.spherical.to_direction()
    }

    pub fn to_earth_equatorial(&self) -> crate::earth_equatorial::EarthEquatorial {
        self.to_direction().to_earth_equatorial()
    }

    pub fn to_equatorial(&self, axis: Direction) -> Equatorial {
        self.to_direction().to_equatorial(axis)
    }

    pub fn to_spherical(&self) -> Spherical {
        self.spherical
    }

    pub fn angle_to(&self, other: &Self) -> Angle {
        self.to_direction().angle_to(&other.to_direction())
    }
}

impl Neg for &Ecliptic {
    type Output = Ecliptic;

    fn neg(self) -> Ecliptic {
        Ecliptic {
            spherical: -self.spherical,
        }
    }
}

impl Neg for Ecliptic {
    type Output = Ecliptic;

    fn neg(self) -> Ecliptic {
        -&self
    }
}

impl Display for Ecliptic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.spherical)
    }
}

#[cfg(test)]
pub(super) const EARTH_NORTH_POLE_IN_ECLIPTIC_COORDINATES: Ecliptic =
    Ecliptic::new(Spherical::new(
        crate::angle_helper::QUARTER_CIRC,
        Angle {
            rad: crate::angle_helper::QUARTER_CIRC.get::<radian>()
                - crate::angle_helper::EARTH_AXIS_TILT.get::<radian>(),
        },
    ));

#[cfg(test)]
mod tests {

    use uom::si::{angle::radian, f64::Angle};

    use crate::{angle_helper::angle_eq_within, direction::Direction, traits::*};

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
                        let angle = Angle::new::<radian>((*angle).abs());
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
