use std::fmt::Display;

use crate::astro_display::AstroDisplay;

use super::{direction::Direction, spherical::SphericalCoordinates};
use simple_si_units::geometry::Angle;

pub struct EquatorialCoordinates {
    spherical: SphericalCoordinates,
    rotation_axis: Direction,
}

impl EquatorialCoordinates {
    pub const fn new(spherical: SphericalCoordinates, rotation_axis: Direction) -> Self {
        Self {
            spherical,
            rotation_axis,
        }
    }

    pub fn get_longitude(&self) -> Angle<f64> {
        self.spherical.get_longitude()
    }

    pub fn get_latitude(&self) -> Angle<f64> {
        self.spherical.get_latitude()
    }

    pub fn set_longitude(&mut self, longitude: Angle<f64>) {
        self.spherical.set_longitude(longitude);
    }

    pub fn set_latitude(&mut self, latitude: Angle<f64>) {
        self.spherical.set_latitude(latitude);
    }

    pub fn get_spherical(&self) -> &SphericalCoordinates {
        &self.spherical
    }

    pub fn get_rotation_axis(&self) -> &Direction {
        &self.rotation_axis
    }

    pub(crate) fn to_direction(&self) -> Direction {
        self.spherical
            .to_direction()
            .active_rotation_to_new_z_axis(&self.rotation_axis)
    }
}

impl AstroDisplay for EquatorialCoordinates {
    fn astro_display(&self) -> String {
        format!(
            "{} relative to equatorial plane with north pole at {}",
            self.spherical, self.rotation_axis
        )
    }
}

impl Display for EquatorialCoordinates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.astro_display())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        coordinates::{
            direction::Direction,
            earth_equatorial::{
                EarthEquatorialCoordinates, EARTH_NORTH_POLE_IN_ECLIPTIC_COORDINATES,
            },
            equatorial::EquatorialCoordinates,
            spherical::SphericalCoordinates,
        },
        tests::TEST_ACCURACY,
    };

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

                    let coordinates =
                        EquatorialCoordinates::new(SphericalCoordinates::Z_DIRECTION, axis.clone());
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

                    let coordinates = EquatorialCoordinates::new(
                        -SphericalCoordinates::Z_DIRECTION,
                        axis.clone(),
                    );
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

                    let coordinates =
                        EquatorialCoordinates::new(SphericalCoordinates::X_DIRECTION, axis);
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

                    let coordinates =
                        EquatorialCoordinates::new(-SphericalCoordinates::X_DIRECTION, axis);
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
            .get_spherical()
            .to_direction();

        for long in ordinates.clone() {
            for lat in ordinates.clone() {
                let long = Angle::from_radians(long);
                let lat = Angle::from_radians(lat);
                let spherical = SphericalCoordinates::new(long, lat);

                let equatorial_coordinates =
                    EquatorialCoordinates::new(spherical, earth_north.clone());
                let earth_equatorial_coordinates = EarthEquatorialCoordinates::new(long, lat);

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

        let equatorial_x =
            EquatorialCoordinates::new(SphericalCoordinates::X_DIRECTION, axis.clone());
        let expected = Direction::X;
        let actual = equatorial_x.to_direction();
        println!("expected: {},\n actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));

        let equatorial_y =
            EquatorialCoordinates::new(SphericalCoordinates::Y_DIRECTION, axis.clone());
        let expected = -&Direction::Z;
        let actual = equatorial_y.to_direction();
        println!("expected: {},\n actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));
    }

    #[test]
    fn axis_tilted_to_x() {
        let axis = Direction::X;

        let equatorial_x =
            EquatorialCoordinates::new(SphericalCoordinates::X_DIRECTION, axis.clone());
        let expected = -&Direction::Y;
        let actual = equatorial_x.to_direction();
        println!("expected: {},\n actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));

        let equatorial_y =
            EquatorialCoordinates::new(SphericalCoordinates::Y_DIRECTION, axis.clone());
        let expected = -&Direction::Z;
        let actual = equatorial_y.to_direction();
        println!("expected: {},\n actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, TEST_ACCURACY));
    }
}
