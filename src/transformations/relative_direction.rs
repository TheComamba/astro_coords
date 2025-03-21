//! Functions for transforming directions relative to a normal vector.

use crate::direction::Direction;
use crate::traits::*;

/// Rotates object_direction as it would seem to an observer looking in the direction of observer_normal.
///
/// The rotation_reference is used to determine the direction of the rotation.
/// Its projection onto the new XY plane points along the new X axis.
pub fn direction_relative_to_normal(
    object_direction: &Direction,
    observer_normal: &Direction,
    rotation_reference: &Direction,
) -> Direction {
    let object_direction = object_direction.passive_rotation_to_new_z_axis(observer_normal);
    let rotation_reference = rotation_reference.passive_rotation_to_new_z_axis(observer_normal);
    let new_x = Direction::new(rotation_reference.x(), rotation_reference.y(), 0.);
    let new_x = match new_x {
        Ok(new_x) => new_x,
        Err(_) => return object_direction,
    };
    let mut rotation_angle = new_x.angle_to(&Direction::X);
    if new_x.y() > 0. {
        rotation_angle = -rotation_angle;
    }
    object_direction.rotated(rotation_angle, &Direction::Z)
}

#[cfg(test)]
mod tests {
    use uom::si::{angle::radian, f64::Angle};

    use crate::angle_helper::*;

    use super::*;

    const TEST_ACCURACY: f64 = 1e-5;

    #[test]
    fn direction_relative_to_itself_is_z_axis() {
        let ordinates = vec![0., 1., -1., 12.];
        for x in ordinates.clone().iter() {
            for y in ordinates.clone().iter() {
                for z in ordinates.clone().iter() {
                    let dir = Direction::new(*x, *y, *z);
                    if dir.is_err() {
                        continue;
                    }
                    let dir = dir.unwrap();
                    let rotation_reference = dir.some_orthogonal_vector();
                    let rotated = direction_relative_to_normal(&dir, &dir, &rotation_reference);

                    assert!(rotated.eq_within(&Direction::Z, TEST_ACCURACY));
                }
            }
        }
    }

    #[test]
    fn objects_along_the_new_x_axis_have_no_y_component() {
        let ordinates = vec![0., 1., -1., 12.];
        let angles = vec![
            angle_zero(),
            Angle::new::<radian>(1.),
            Angle::new::<radian>(-1.),
            Angle::new::<radian>(12.),
        ];
        for x1 in ordinates.clone().iter() {
            for y1 in ordinates.clone().iter() {
                for z1 in ordinates.clone().iter() {
                    for x2 in ordinates.clone().iter() {
                        for y2 in ordinates.clone().iter() {
                            for z2 in ordinates.clone().iter() {
                                for angle in angles.iter() {
                                    let observer_normal = Direction::new(*x1, *y1, *z1);
                                    let rotation_reference = Direction::new(*x2, *y2, *z2);
                                    if observer_normal.is_err() || rotation_reference.is_err() {
                                        continue;
                                    }
                                    let observer_normal = observer_normal.unwrap();
                                    let rotation_reference = rotation_reference.unwrap();
                                    let ortho = observer_normal.cross_product(&rotation_reference);
                                    if ortho.is_err() {
                                        continue;
                                    }
                                    let ortho = ortho.unwrap();
                                    let object_direction =
                                        rotation_reference.rotated(*angle, &ortho);

                                    let dir_in_new_system = direction_relative_to_normal(
                                        &object_direction,
                                        &observer_normal,
                                        &rotation_reference,
                                    );

                                    println!("\nobject_direction: {}", object_direction);
                                    println!("observer_normal: {}", observer_normal);
                                    println!("rotation_reference: {}", rotation_reference);
                                    println!("dir_in_new_system: {}", dir_in_new_system);
                                    println!("y component: {}", dir_in_new_system.y());
                                    assert!(dir_in_new_system.y().abs() < 10. * TEST_ACCURACY);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn rotating_to_original_x_and_z_does_nothing() {
        let ordinates = vec![0., 1., -1., 12.];
        for x in ordinates.clone().iter() {
            for y in ordinates.clone().iter() {
                for z in ordinates.clone().iter() {
                    let object_direction = Direction::new(*x, *y, *z);
                    if object_direction.is_err() {
                        continue;
                    }
                    let orig = object_direction.unwrap();

                    let new = direction_relative_to_normal(&orig, &Direction::Z, &Direction::X);

                    assert!(new.eq_within(&orig, TEST_ACCURACY));
                }
            }
        }
    }

    #[test]
    fn mirroring_x_axis_rotates_around_z() {
        let ordinates = vec![0., 1., -1., 12.];
        for x in ordinates.clone().iter() {
            for y in ordinates.clone().iter() {
                for z in ordinates.clone().iter() {
                    let object_direction = Direction::new(*x, *y, *z);
                    if object_direction.is_err() {
                        continue;
                    }
                    let orig = object_direction.unwrap();

                    let new = direction_relative_to_normal(&orig, &Direction::Z, &-&Direction::X);

                    let expected = orig.rotated(half_circ(), &Direction::Z);

                    assert!(new.eq_within(&expected, TEST_ACCURACY));
                }
            }
        }
    }

    #[test]
    fn mirroring_z_axis_rotates_around_x() {
        let ordinates = vec![0., 1., -1., 12.];
        for x in ordinates.clone().iter() {
            for y in ordinates.clone().iter() {
                for z in ordinates.clone().iter() {
                    let object_direction = Direction::new(*x, *y, *z);
                    if object_direction.is_err() {
                        continue;
                    }
                    let orig = object_direction.unwrap();

                    let new = direction_relative_to_normal(&orig, &-&Direction::Z, &Direction::X);

                    let expected = orig.rotated(half_circ(), &Direction::X);

                    assert!(new.eq_within(&expected, TEST_ACCURACY));
                }
            }
        }
    }

    #[test]
    fn converting_axes() {
        let new = direction_relative_to_normal(&Direction::X, &Direction::Y, &Direction::Z);
        let expected = Direction::Y;
        assert!(new.eq_within(&expected, TEST_ACCURACY));

        let new = direction_relative_to_normal(&Direction::X, &Direction::Z, &Direction::Y);
        let expected = -&Direction::Y;
        assert!(new.eq_within(&expected, TEST_ACCURACY));

        let new = direction_relative_to_normal(&Direction::Y, &Direction::X, &Direction::Z);
        let expected = -&Direction::Y;
        assert!(new.eq_within(&expected, TEST_ACCURACY));

        let new = direction_relative_to_normal(&Direction::Y, &Direction::Z, &Direction::X);
        let expected = Direction::Y;
        assert!(new.eq_within(&expected, TEST_ACCURACY));

        let new = direction_relative_to_normal(&Direction::Z, &Direction::X, &Direction::Y);
        let expected = Direction::Y;
        assert!(new.eq_within(&expected, TEST_ACCURACY));

        let new = direction_relative_to_normal(&Direction::Z, &Direction::Y, &Direction::X);
        let expected = -&Direction::Y;
        assert!(new.eq_within(&expected, TEST_ACCURACY));
    }
}
