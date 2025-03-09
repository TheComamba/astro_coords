//! This module contains the Equatorial struct and its implementation.

use std::fmt::Display;

use uom::si::f64::{Angle, Length};

use crate::{
    cartesian::Cartesian, earth_equatorial::EarthEquatorial, ecliptic::Ecliptic, traits::*,
};

use super::{direction::Direction, spherical::Spherical};

pub struct Equatorial {
    pub spherical: Spherical,
    pub rotation_axis: Direction,
}

impl Equatorial {
    pub const fn new(spherical: Spherical, rotation_axis: Direction) -> Self {
        Self {
            spherical,
            rotation_axis,
        }
    }

    pub fn to_cartesian(&self, length: Length) -> Cartesian {
        self.to_direction().to_cartesian(length)
    }

    pub fn to_direction(&self) -> Direction {
        self.spherical
            .to_direction()
            .active_rotation_to_new_z_axis(&self.rotation_axis)
    }

    pub fn to_earth_equatorial(&self) -> EarthEquatorial {
        self.to_direction().to_earth_equatorial()
    }

    pub fn to_ecliptic(&self) -> Ecliptic {
        self.to_direction().to_ecliptic()
    }

    pub fn to_spherical(&self) -> Spherical {
        self.to_ecliptic().to_spherical()
    }

    pub fn eq_within(&self, other: &Equatorial, accuracy: Angle) -> bool {
        self.to_ecliptic().eq_within(&other.to_ecliptic(), accuracy)
    }
}

impl Display for Equatorial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} relative to equatorial plane with north pole at {}",
            self.spherical, self.rotation_axis
        )
    }
}

#[cfg(test)]
mod tests {
    use uom::si::angle::radian;

    use crate::{
        earth_equatorial::EarthEquatorial, ecliptic::EARTH_NORTH_POLE_IN_ECLIPTIC_COORDINATES,
    };

    use super::*;

    const TEST_ACCURACY: f64 = 1e-5;

    #[test]
    fn north_pole_points_along_axis() {
        let ordinates: Vec<f64> = vec![-1., 0., 1., 10.];
        for x in ordinates.clone() {
            for y in ordinates.clone() {
                for z in ordinates.clone() {
                    let axis = Direction::new(x, y, z);
                    if axis.is_err() {
                        continue;
                    }
                    let axis = axis.unwrap();

                    let coordinates = Equatorial::new(Spherical::Z_DIRECTION, axis.clone());
                    let expected = axis;
                    let actual = coordinates.to_direction();
                    println!("expected: {},\n actual: {}", expected, actual);
                    assert!(actual.eq_within(&expected, TEST_ACCURACY));
                }
            }
        }
    }

    #[test]
    fn south_pole_points_along_negative_axis() {
        let ordinates: Vec<f64> = vec![-1., 0., 1., 10.];
        for x in ordinates.clone() {
            for y in ordinates.clone() {
                for z in ordinates.clone() {
                    let axis = Direction::new(x, y, z);
                    if axis.is_err() {
                        continue;
                    }
                    let axis = axis.unwrap();

                    let coordinates = Equatorial::new(-Spherical::Z_DIRECTION, axis.clone());
                    let expected = -&axis;
                    let actual = coordinates.to_direction();
                    println!("expected: {},\n actual: {}", expected, actual);
                    assert!(actual.eq_within(&expected, TEST_ACCURACY));
                }
            }
        }
    }

    #[test]
    fn x_axis_lies_in_horizontal_plane() {
        let ordinates: Vec<f64> = vec![-1., 0., 1., 10.];
        for x in ordinates.clone() {
            for y in ordinates.clone() {
                for z in ordinates.clone() {
                    let axis = Direction::new(x, y, z);
                    if axis.is_err() {
                        continue;
                    }
                    let axis = axis.unwrap();

                    let coordinates = Equatorial::new(Spherical::X_DIRECTION, axis);
                    let direction = coordinates.to_direction();
                    assert!(direction.z().abs() < TEST_ACCURACY);
                }
            }
        }
    }

    #[test]
    fn minus_x_axis_lies_in_horizontal_plane() {
        let ordinates: Vec<f64> = vec![-1., 0., 1., 10.];
        for x in ordinates.clone() {
            for y in ordinates.clone() {
                for z in ordinates.clone() {
                    let axis = Direction::new(x, y, z);
                    if axis.is_err() {
                        continue;
                    }
                    let axis = axis.unwrap();

                    let coordinates = Equatorial::new(-Spherical::X_DIRECTION, axis);
                    let direction = coordinates.to_direction();
                    assert!(direction.z().abs() < TEST_ACCURACY);
                }
            }
        }
    }

    #[test]
    fn behaves_like_earth_equatorial() {
        let ordinates: Vec<f64> = vec![-1., 0., 1., 10.];
        let earth_north = EARTH_NORTH_POLE_IN_ECLIPTIC_COORDINATES
            .spherical
            .to_direction();

        for long in ordinates.clone() {
            for lat in ordinates.clone() {
                let long = Angle::new::<radian>(long);
                let lat = Angle::new::<radian>(lat);
                let spherical = Spherical::new(long, lat);

                let equatorial_coordinates = Equatorial::new(spherical, earth_north.clone());
                let earth_equatorial_coordinates = EarthEquatorial::new(long, lat);

                let expected = earth_equatorial_coordinates.to_direction();
                let actual = equatorial_coordinates.to_direction();
                println!("expected: {},\n actual: {}", expected, actual);
                assert!(actual.eq_within(&expected, TEST_ACCURACY));
            }
        }
    }

    #[test]
    fn axis_tilted_to_y() {
        let axis = Direction::Y;

        let equatorial_x = Equatorial::new(Spherical::X_DIRECTION, axis.clone());
        let expected = Direction::X;
        let actual = equatorial_x.to_direction();
        println!("expected: {},\n actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));

        let equatorial_y = Equatorial::new(Spherical::Y_DIRECTION, axis.clone());
        let expected = -&Direction::Z;
        let actual = equatorial_y.to_direction();
        println!("expected: {},\n actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));
    }

    #[test]
    fn axis_tilted_to_x() {
        let axis = Direction::X;

        let equatorial_x = Equatorial::new(Spherical::X_DIRECTION, axis.clone());
        let expected = -&Direction::Y;
        let actual = equatorial_x.to_direction();
        println!("expected: {},\n actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));

        let equatorial_y = Equatorial::new(Spherical::Y_DIRECTION, axis.clone());
        let expected = -&Direction::Z;
        let actual = equatorial_y.to_direction();
        println!("expected: {},\n actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));
    }
}
