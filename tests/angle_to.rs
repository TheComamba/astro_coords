use uom::si::angle::radian;
use utils::{constants::*, examples::*};

mod utils;

#[test]
fn angle_is_the_same_for_direction_and_cartesian() {
    for dir1 in direction_examples() {
        for dir2 in direction_examples() {
            for dist1 in positive_distance_examples() {
                for dist2 in positive_distance_examples() {
                    let angle = dir1.angle_to(&dir2);
                    let cart1 = dir1.to_cartesian(dist1);
                    let cart2 = dir2.to_cartesian(dist2);
                    let new_angle = cart1.angle_to(&cart2).unwrap();
                    assert!((angle - new_angle).get::<radian>().abs() < ACC,
                    "dir1: {:?}\ndir2: {:?}\ndist1: {:?}\ndist2: {:?}\nangle: {:?}\nnew_angle: {:?}",
                    dir1, dir2, dist1, dist2, angle, new_angle);
                }
            }
        }
    }
}

#[test]
fn angle_is_the_same_for_direction_and_spherical() {
    for dir1 in direction_examples() {
        for dir2 in direction_examples() {
            let angle = dir1.angle_to(&dir2);
            let sph1 = dir1.to_spherical();
            let sph2 = dir2.to_spherical();
            let new_angle = sph1.angle_to(&sph2);
            assert!(
                (angle - new_angle).get::<radian>().abs() < ACC,
                "dir1: {:?}\ndir2: {:?}\nsph1: {:?}\nsph2: {:?}\nangle: {:?}\nnew_angle: {:?}",
                dir1,
                dir2,
                sph1,
                sph2,
                angle,
                new_angle
            );
        }
    }
}

#[test]
fn angle_is_the_same_for_direction_and_spherical2() {
    for sph1 in spherical_examples() {
        for sph2 in spherical_examples() {
            let angle = sph1.angle_to(&sph2);
            let dir1 = sph1.to_direction();
            let dir2 = sph2.to_direction();
            let new_angle = dir1.angle_to(&dir2);
            assert!(
                (angle - new_angle).get::<radian>().abs() < ACC,
                "dir1: {:?}\ndir2: {:?}\nsph1: {:?}\nsph2: {:?}\nangle: {:?}\nnew_angle: {:?}",
                dir1,
                dir2,
                sph1,
                sph2,
                angle,
                new_angle
            );
        }
    }
}
