use astro_coords::{cartesian::Cartesian, direction::Direction, spherical::Spherical, traits::*};
use simple_si_units::geometry::Angle;
use utils::{constants::*, examples::*};

mod utils;

fn cartesian_points_along_z(cart: &Cartesian) -> bool {
    cart.x.m.abs() < ACC && cart.y.m.abs() < ACC
}

#[test]
fn cartesian_roundtrips() {
    for target_frame in reference_frame_examples() {
        for physical in physical_cartesian_examples() {
            let original_frame = physical.reference_frame();
            let original_mathematical = physical.mathematical_coordinates();
            let mut new_physical = physical.clone();

            new_physical.change_reference_frame(target_frame);
            assert_eq!(new_physical.reference_frame(), target_frame);
            new_physical.change_reference_frame(original_frame);

            let mathematical_after_roundtrip = new_physical.mathematical_coordinates();
            assert!(
                original_mathematical.eq_within(mathematical_after_roundtrip, DISTANCE_ACC),
                "{:?} != {:?}",
                original_mathematical,
                mathematical_after_roundtrip
            );
        }
    }
}

fn direction_points_along_z(dir: &Direction) -> bool {
    dir.x().abs() < ACC && dir.y().abs() < ACC
}

#[test]
fn direction_roundtrips() {
    for target_frame in reference_frame_examples() {
        for physical in physical_direction_examples() {
            let original_frame = physical.reference_frame();
            let original_mathematical = physical.mathematical_coordinates();
            let mut new_physical = physical.clone();

            new_physical.change_reference_frame(target_frame);
            assert_eq!(new_physical.reference_frame(), target_frame);
            new_physical.change_reference_frame(original_frame);

            let mathematical_after_roundtrip = new_physical.mathematical_coordinates();
            assert!(
                original_mathematical.eq_within(mathematical_after_roundtrip, ACC),
                "{:?} != {:?}",
                original_mathematical,
                mathematical_after_roundtrip
            );
        }
    }
}

fn spherical_points_along_z(sph: &Spherical) -> bool {
    (sph.latitude - Angle::from_deg(90.)).rad.abs() < ACC
}

#[test]
fn spherical_roundtrips() {
    for target_frame in reference_frame_examples() {
        for physical in physical_spherical_examples() {
            let original_frame = physical.reference_frame();
            let original_mathematical = physical.mathematical_coordinates();
            let mut new_physical = physical.clone();

            new_physical.change_reference_frame(target_frame);
            assert_eq!(new_physical.reference_frame(), target_frame);
            new_physical.change_reference_frame(original_frame);

            let mathematical_after_roundtrip = new_physical.mathematical_coordinates();
            assert!(
                original_mathematical.eq_within(mathematical_after_roundtrip, ANGLE_ACC),
                "{:?} != {:?}",
                original_mathematical,
                mathematical_after_roundtrip
            );
        }
    }
}

#[test]
fn changing_reference_frame_of_cartesian_preserves_length() {
    for target_frame in reference_frame_examples() {
        for physical in physical_cartesian_examples() {
            let original_length = physical.mathematical_coordinates().length();
            let mut new_physical = physical.clone();

            new_physical.change_reference_frame(target_frame);
            let new_length = new_physical.mathematical_coordinates().length();
            assert!(
                (original_length - new_length).m.abs() < ACC,
                "{} != {}",
                original_length,
                new_length
            );
        }
    }
}
