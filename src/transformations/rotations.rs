//! Functions for rotating vectors.

use simple_si_units::geometry::Angle;
use std::ops::{Add, Mul};

use crate::direction::Direction;

pub(crate) fn rotated_tuple<T>(tup: (T, T, T), angle: Angle<f64>, axis: &Direction) -> (T, T, T)
where
    T: Mul<f64, Output = T> + Add<Output = T> + Copy,
{
    let cos = angle.rad.cos();
    let sin = angle.rad.sin();

    let (x, y, z) = tup;

    let ux = axis.x();
    let uy = axis.y();
    let uz = axis.z();

    let r_11 = cos + ux * ux * (1. - cos);
    let r_12 = ux * uy * (1. - cos) - uz * sin;
    let r_13 = ux * uz * (1. - cos) + uy * sin;

    let r_21 = uy * ux * (1. - cos) + uz * sin;
    let r_22 = cos + uy * uy * (1. - cos);
    let r_23 = uy * uz * (1. - cos) - ux * sin;

    let r_31 = uz * ux * (1. - cos) - uy * sin;
    let r_32 = uz * uy * (1. - cos) + ux * sin;
    let r_33 = cos + uz * uz * (1. - cos);

    let x_out = x * r_11 + y * r_12 + z * r_13;
    let y_out = x * r_21 + y * r_22 + z * r_23;
    let z_out = x * r_31 + y * r_32 + z * r_33;
    (x_out, y_out, z_out)
}

/// Calculates the angle and axis of rotation to rotate from one direction to another.
///
/// # Example
/// ```
/// use simple_si_units::geometry::Angle;
/// use astro_coordinates::direction::Direction;
///
/// let start = Direction::new(1., 0., 0.).unwrap();
/// let end = Direction::new(0., 1., 0.).unwrap();
/// let (angle, axis) = astro_coordinates::transformations::get_rotation_parameters(&start, &end);
/// assert_eq!(angle, Angle::from_degrees(90.));
/// assert_eq!(axis, Direction::Z);
/// ```
pub fn get_rotation_parameters(start: &Direction, end: &Direction) -> (Angle<f64>, Direction) {
    let angle = start.angle_to(end);
    let axis = start.cross_product(end);
    if let Ok(axis) = axis {
        (angle, axis)
    } else {
        // start and and end are (anti) parallel
        let axis = start.some_orthogonal_vector();
        (angle, axis)
    }
}

#[cfg(test)]
mod tests {
    use crate::angle_helper::{test::*, *};

    use super::*;

    const TEST_ACCURACY: f64 = 1e-5;

    const X_VECTOR: (f64, f64, f64) = (1., 0., 0.);
    const MINUS_X_VECTOR: (f64, f64, f64) = (-1., 0., 0.);
    const Y_VECTOR: (f64, f64, f64) = (0., 1., 0.);
    const MINUS_Y_VECTOR: (f64, f64, f64) = (0., -1., 0.);
    const Z_VECTOR: (f64, f64, f64) = (0., 0., 1.);
    const MINUS_Z_VECTOR: (f64, f64, f64) = (0., 0., -1.);

    fn kinda_equal(a: (f64, f64, f64), b: (f64, f64, f64)) -> bool {
        let (ax, ay, az) = a;
        let (bx, by, bz) = b;
        (ax - bx).abs() < TEST_ACCURACY
            && (ay - by).abs() < TEST_ACCURACY
            && (az - bz).abs() < TEST_ACCURACY
    }

    fn print_expectations(expected: (f64, f64, f64), actual: (f64, f64, f64)) {
        println!(
            "expected ({}, {}, {}), actual ({}, {}, {})",
            expected.0, expected.1, expected.2, actual.0, actual.1, actual.2
        );
    }

    #[test]
    fn test_rotating_x_around_z() {
        let start = X_VECTOR;

        let rotated = rotated_tuple(start, QUARTER_CIRC, &Direction::Z);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, HALF_CIRC, &Direction::Z);
        let expected = MINUS_X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, THREE_QUARTER_CIRC, &Direction::Z);
        let expected = MINUS_Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, FULL_CIRC, &Direction::Z);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_y_around_z() {
        let start = Y_VECTOR;

        let rotated = rotated_tuple(start, QUARTER_CIRC, &Direction::Z);
        let expected = MINUS_X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, HALF_CIRC, &Direction::Z);
        let expected = MINUS_Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, THREE_QUARTER_CIRC, &Direction::Z);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, FULL_CIRC, &Direction::Z);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_z_around_z() {
        let start = Z_VECTOR;

        let rotated = rotated_tuple(start, QUARTER_CIRC, &Direction::Z);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, HALF_CIRC, &Direction::Z);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, THREE_QUARTER_CIRC, &Direction::Z);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, FULL_CIRC, &Direction::Z);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_x_around_x() {
        let start = X_VECTOR;

        let rotated = rotated_tuple(start, QUARTER_CIRC, &Direction::X);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, HALF_CIRC, &Direction::X);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, THREE_QUARTER_CIRC, &Direction::X);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, FULL_CIRC, &Direction::X);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_y_around_x() {
        let start = Y_VECTOR;

        let rotated = rotated_tuple(start, QUARTER_CIRC, &Direction::X);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, HALF_CIRC, &Direction::X);
        let expected = MINUS_Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, THREE_QUARTER_CIRC, &Direction::X);
        let expected = MINUS_Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, FULL_CIRC, &Direction::X);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_z_around_x() {
        let start = Z_VECTOR;

        let rotated = rotated_tuple(start, QUARTER_CIRC, &Direction::X);
        let expected = MINUS_Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, HALF_CIRC, &Direction::X);
        let expected = MINUS_Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, THREE_QUARTER_CIRC, &Direction::X);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, FULL_CIRC, &Direction::X);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_x_around_y() {
        let start = X_VECTOR;

        let rotated = rotated_tuple(start, QUARTER_CIRC, &Direction::Y);
        let expected = MINUS_Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, HALF_CIRC, &Direction::Y);
        let expected = MINUS_X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, THREE_QUARTER_CIRC, &Direction::Y);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, FULL_CIRC, &Direction::Y);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_y_around_y() {
        let start = Y_VECTOR;

        let rotated = rotated_tuple(start, QUARTER_CIRC, &Direction::Y);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, HALF_CIRC, &Direction::Y);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, THREE_QUARTER_CIRC, &Direction::Y);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, FULL_CIRC, &Direction::Y);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_z_around_y() {
        let start = Z_VECTOR;

        let rotated = rotated_tuple(start, QUARTER_CIRC, &Direction::Y);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, HALF_CIRC, &Direction::Y);
        let expected = MINUS_Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, THREE_QUARTER_CIRC, &Direction::Y);
        let expected = MINUS_X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, FULL_CIRC, &Direction::Y);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_around_diagonal_axis() {
        let axis = Direction::new(1., 1., 1.).unwrap();

        let rotated = rotated_tuple(X_VECTOR, ONE_THIRD_CIRC, &axis);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(Y_VECTOR, ONE_THIRD_CIRC, &axis);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(Z_VECTOR, ONE_THIRD_CIRC, &axis);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(X_VECTOR, TWO_THIRDS_CIRC, &axis);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(Y_VECTOR, TWO_THIRDS_CIRC, &axis);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(Z_VECTOR, TWO_THIRDS_CIRC, &axis);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(X_VECTOR, -ONE_THIRD_CIRC, &axis);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(Y_VECTOR, -ONE_THIRD_CIRC, &axis);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(Z_VECTOR, -ONE_THIRD_CIRC, &axis);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn get_rotation_parameters_test() {
        const ROTATION_DIRECTION_ACCURACY: f64 = 1e-3;
        let ordinates = vec![-1., 0., 1., 10.];
        for start_x in ordinates.clone().iter() {
            for start_y in ordinates.clone().iter() {
                for start_z in ordinates.clone().iter() {
                    for end_x in ordinates.clone().iter() {
                        for end_y in ordinates.clone().iter() {
                            for end_z in ordinates.clone().iter() {
                                let start = Direction::new(*start_x, *start_y, *start_z);
                                let end = Direction::new(*end_x, *end_y, *end_z);
                                if start.is_err() || end.is_err() {
                                    continue;
                                }
                                let start = start.unwrap();
                                let end = end.unwrap();
                                println!("start: {}, end: {}", start, end);
                                let (angle, axis) = get_rotation_parameters(&start, &end);
                                println!("angle: {}, axis: {}", angle, axis);
                                let rotated = start.rotated(angle, &axis);
                                println!("expected: {}, actual: {}", end, rotated);
                                assert!(rotated.eq_within(&end, ROTATION_DIRECTION_ACCURACY));
                            }
                        }
                    }
                }
            }
        }
    }
}
