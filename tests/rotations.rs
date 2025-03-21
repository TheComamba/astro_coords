use astro_coords::{direction::Direction, traits::*};
use uom::{fmt::DisplayStyle, si::angle::radian};
use utils::{constants::*, examples::*};

mod utils;

#[test]
fn rotated_x_is_the_same_as_rotation_around_x() {
    for vec in direction_examples() {
        for angle in angle_examples() {
            let rotated_1 = vec.rotated_x(angle);
            let rotated_2 = vec.rotated(angle, &Direction::X);
            assert!(rotated_1.eq_within(&rotated_2, ACC));
        }
    }
}

#[test]
fn rotated_y_is_the_same_as_rotation_around_y() {
    for vec in direction_examples() {
        for angle in angle_examples() {
            let rotated_1 = vec.rotated_y(angle);
            let rotated_2 = vec.rotated(angle, &Direction::Y);
            assert!(rotated_1.eq_within(&rotated_2, ACC));
        }
    }
}

#[test]
fn rotated_z_is_the_same_as_rotation_around_z() {
    for vec in direction_examples() {
        for angle in angle_examples() {
            let rotated_1 = vec.rotated_z(angle);
            let rotated_2 = vec.rotated(angle, &Direction::Z);
            assert!(rotated_1.eq_within(&rotated_2, ACC));
        }
    }
}

#[test]
fn rotated_is_the_same_for_direction_and_spherical() {
    for vec in direction_examples() {
        for angle in angle_examples() {
            for axis in direction_examples() {
                let rotated_dir = vec.rotated(angle, &axis);

                let spherical_dir = vec.to_spherical();
                let rotated_spherical = spherical_dir.rotated(angle, &axis);

                let expected = rotated_dir.to_spherical();
                assert!(
                    rotated_spherical.eq_within(&expected, ANGLE_ACC()),
                    "{} != {}",
                    rotated_spherical,
                    expected
                );
            }
        }
    }
}

#[test]
fn rotated_x_is_the_same_for_direction_and_spherical() {
    for vec in direction_examples() {
        for angle in angle_examples() {
            let rotated_dir = vec.rotated_x(angle);

            let spherical_dir = vec.to_spherical();
            let rotated_spherical = spherical_dir.rotated_x(angle);

            let expected = rotated_dir.to_spherical();
            assert!(
                rotated_spherical.eq_within(&expected, ANGLE_ACC()),
                "{} != {}",
                rotated_spherical,
                expected
            );
        }
    }
}

#[test]
fn rotated_y_is_the_same_for_direction_and_spherical() {
    for vec in direction_examples() {
        for angle in angle_examples() {
            let rotated_dir = vec.rotated_y(angle);

            let spherical_dir = vec.to_spherical();
            let rotated_spherical = spherical_dir.rotated_y(angle);

            let expected = rotated_dir.to_spherical();
            assert!(
                rotated_spherical.eq_within(&expected, ANGLE_ACC()),
                "{} != {}",
                rotated_spherical,
                expected
            );
        }
    }
}

#[test]
fn rotated_z_is_the_same_for_direction_and_spherical() {
    for vec in direction_examples() {
        for angle in angle_examples() {
            let rotated_dir = vec.rotated_z(angle);

            let spherical_dir = vec.to_spherical();
            let rotated_spherical = spherical_dir.rotated_z(angle);

            let expected = rotated_dir.to_spherical();
            assert!(
                rotated_spherical.eq_within(&expected, ANGLE_ACC()),
                "vec: {}\nangle: {}\n{} != {}",
                vec,
                angle.into_format_args(radian, DisplayStyle::Abbreviation),
                rotated_spherical,
                expected
            );
        }
    }
}

#[test]
fn active_rotation_is_the_same_for_direction_and_spherical() {
    for new_z in direction_examples() {
        for vec in direction_examples() {
            let rotated_dir = vec.active_rotation_to_new_z_axis(&new_z);

            let spherical_new_z = new_z.to_spherical();
            let spherical_dir = vec.to_spherical();
            let rotated_spherical = spherical_dir.active_rotation_to_new_z_axis(&spherical_new_z);

            let expected = rotated_dir.to_spherical();
            assert!(
                rotated_spherical.eq_within(&expected, ANGLE_ACC()),
                "new_z: {}\nvec: {}\nrotated_dir: {}\nrotated_spherical: {}\nexpected: {}",
                new_z,
                vec,
                rotated_dir,
                rotated_spherical,
                expected
            );
        }
    }
}

#[test]
fn passive_rotation_is_the_same_for_direction_and_spherical() {
    for new_z in direction_examples() {
        for vec in direction_examples() {
            let rotated_dir = vec.passive_rotation_to_new_z_axis(&new_z);

            let spherical_new_z = new_z.to_spherical();
            let spherical_dir = vec.to_spherical();
            let rotated_spherical = spherical_dir.passive_rotation_to_new_z_axis(&spherical_new_z);

            let expected = rotated_dir.to_spherical();
            assert!(rotated_spherical.eq_within(&expected, ANGLE_ACC()));
        }
    }
}
