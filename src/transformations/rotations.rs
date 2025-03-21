//! Functions for rotating vectors.

use std::ops::{Add, Mul, Sub};

use uom::si::{angle::radian, f64::Angle};

use crate::direction::Direction;

/// Rotates a 3-tuple around an arbitrary axis.
pub(crate) fn rotated_tuple<T>(tup: (T, T, T), angle: Angle, axis: &Direction) -> (T, T, T)
where
    T: Mul<f64, Output = T> + Add<Output = T> + Copy,
{
    let cos = angle.get::<radian>().cos();
    let sin = angle.get::<radian>().sin();

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

/// Rotates a 3-tuple around the x-axis.
pub(crate) fn rotated_x_tuple<T>(tup: (T, T, T), angle: Angle) -> (T, T, T)
where
    T: Mul<f64, Output = T> + Add<Output = T> + Sub<Output = T> + Copy,
{
    let cos = angle.get::<radian>().cos();
    let sin = angle.get::<radian>().sin();

    let (x, y, z) = tup;

    let x_out = x;
    let y_out = y * cos - z * sin;
    let z_out = y * sin + z * cos;
    (x_out, y_out, z_out)
}

/// Rotates a 3-tuple around the y-axis.
pub(crate) fn rotated_y_tuple<T>(tup: (T, T, T), angle: Angle) -> (T, T, T)
where
    T: Mul<f64, Output = T> + Add<Output = T> + Sub<Output = T> + Copy,
{
    let cos = angle.get::<radian>().cos();
    let sin = angle.get::<radian>().sin();

    let (x, y, z) = tup;

    let x_out = x * cos + z * sin;
    let y_out = y;
    let z_out = z * cos - x * sin;
    (x_out, y_out, z_out)
}

/// Rotates a 3-tuple around the z-axis.
pub(crate) fn rotated_z_tuple<T>(tup: (T, T, T), angle: Angle) -> (T, T, T)
where
    T: Mul<f64, Output = T> + Add<Output = T> + Sub<Output = T> + Copy,
{
    let cos = angle.get::<radian>().cos();
    let sin = angle.get::<radian>().sin();

    let (x, y, z) = tup;

    let x_out = x * cos - y * sin;
    let y_out = x * sin + y * cos;
    let z_out = z;
    (x_out, y_out, z_out)
}

/// Calculates the angle and axis of rotation to rotate from one direction to another.
///
/// # Example
/// ```
/// use uom::si::f64::Angle;
/// use astro_coords::direction::Direction;
/// use astro_coords::transformations::rotations::get_rotation_parameters;
///
/// let start = Direction::new(1., 0., 0.).unwrap();
/// let end = Direction::new(0., 1., 0.).unwrap();
/// let (angle, axis) = get_rotation_parameters(&start, &end);
/// assert!((angle.get::<degree>() - 90.).abs() < 1e-5);
/// assert!(axis.eq_within(&Direction::Z, 1e-5));
/// ```
pub fn get_rotation_parameters(start: &Direction, end: &Direction) -> (Angle, Direction) {
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
    use crate::traits::*;

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

        let rotated = rotated_tuple(start, quarter_circ(), &Direction::Z);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, half_circ(), &Direction::Z);
        let expected = MINUS_X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, three_quarter_circ(), &Direction::Z);
        let expected = MINUS_Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, full_circ(), &Direction::Z);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_y_around_z() {
        let start = Y_VECTOR;

        let rotated = rotated_tuple(start, quarter_circ(), &Direction::Z);
        let expected = MINUS_X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, half_circ(), &Direction::Z);
        let expected = MINUS_Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, three_quarter_circ(), &Direction::Z);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, full_circ(), &Direction::Z);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_z_around_z() {
        let start = Z_VECTOR;

        let rotated = rotated_tuple(start, quarter_circ(), &Direction::Z);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, half_circ(), &Direction::Z);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, three_quarter_circ(), &Direction::Z);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, full_circ(), &Direction::Z);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_x_around_x() {
        let start = X_VECTOR;

        let rotated = rotated_tuple(start, quarter_circ(), &Direction::X);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, half_circ(), &Direction::X);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, three_quarter_circ(), &Direction::X);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, full_circ(), &Direction::X);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_y_around_x() {
        let start = Y_VECTOR;

        let rotated = rotated_tuple(start, quarter_circ(), &Direction::X);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, half_circ(), &Direction::X);
        let expected = MINUS_Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, three_quarter_circ(), &Direction::X);
        let expected = MINUS_Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, full_circ(), &Direction::X);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_z_around_x() {
        let start = Z_VECTOR;

        let rotated = rotated_tuple(start, quarter_circ(), &Direction::X);
        let expected = MINUS_Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, half_circ(), &Direction::X);
        let expected = MINUS_Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, three_quarter_circ(), &Direction::X);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, full_circ(), &Direction::X);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_x_around_y() {
        let start = X_VECTOR;

        let rotated = rotated_tuple(start, quarter_circ(), &Direction::Y);
        let expected = MINUS_Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, half_circ(), &Direction::Y);
        let expected = MINUS_X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, three_quarter_circ(), &Direction::Y);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, full_circ(), &Direction::Y);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_y_around_y() {
        let start = Y_VECTOR;

        let rotated = rotated_tuple(start, quarter_circ(), &Direction::Y);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, half_circ(), &Direction::Y);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, three_quarter_circ(), &Direction::Y);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, full_circ(), &Direction::Y);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_z_around_y() {
        let start = Z_VECTOR;

        let rotated = rotated_tuple(start, quarter_circ(), &Direction::Y);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, half_circ(), &Direction::Y);
        let expected = MINUS_Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, three_quarter_circ(), &Direction::Y);
        let expected = MINUS_X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(start, full_circ(), &Direction::Y);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));
    }

    #[test]
    fn test_rotating_around_diagonal_axis() {
        let axis = Direction::new(1., 1., 1.).unwrap();

        let rotated = rotated_tuple(X_VECTOR, one_third_circ(), &axis);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(Y_VECTOR, one_third_circ(), &axis);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(Z_VECTOR, one_third_circ(), &axis);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(X_VECTOR, two_thirds_circ(), &axis);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(Y_VECTOR, two_thirds_circ(), &axis);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(Z_VECTOR, two_thirds_circ(), &axis);
        let expected = Y_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(X_VECTOR, -one_third_circ(), &axis);
        let expected = Z_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(Y_VECTOR, -one_third_circ(), &axis);
        let expected = X_VECTOR;
        print_expectations(expected, rotated);
        assert!(kinda_equal(rotated, expected));

        let rotated = rotated_tuple(Z_VECTOR, -one_third_circ(), &axis);
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
                                println!(
                                    "angle: {}, axis: {}",
                                    angle.into_format_args(
                                        radian,
                                        uom::fmt::DisplayStyle::Abbreviation
                                    ),
                                    axis
                                );
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
