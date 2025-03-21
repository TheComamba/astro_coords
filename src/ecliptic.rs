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
    #[inline]
    pub fn x_direction() -> Ecliptic {
        Ecliptic {
            spherical: Spherical::x_direction(),
        }
    }
    #[inline]
    pub fn y_direction() -> Ecliptic {
        Ecliptic {
            spherical: Spherical::y_direction(),
        }
    }
    #[inline]
    pub fn z_direction() -> Ecliptic {
        Ecliptic {
            spherical: Spherical::z_direction(),
        }
    }

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
#[inline]
pub(super) fn earth_north_pole_in_ecliptic_coordinates() -> Ecliptic {
    Ecliptic::new(Spherical::new(
        crate::angle_helper::quarter_circ(),
        crate::angle_helper::quarter_circ() - crate::angle_helper::earth_axis_tilt(),
    ))
}

#[cfg(test)]
mod tests {

    use uom::{
        fmt::DisplayStyle,
        si::{angle::radian, f64::Angle},
    };

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

                        println!(
                            "Expected: {}",
                            angle.into_format_args(radian, DisplayStyle::Abbreviation)
                        );
                        println!(
                            "Actual: {}",
                            actual_angle.into_format_args(radian, DisplayStyle::Abbreviation)
                        );
                        assert!(angle_eq_within(
                            actual_angle,
                            angle,
                            Angle::new::<radian>(1e-5)
                        ));
                    }
                }
            }
        }
    }
}
