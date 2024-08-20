use astro_coords::traits::*;
use utils::{constants::*, examples::*};

mod utils;

#[test]
fn cartesian_roundtrips() {
    for target_frame in reference_frame_examples() {
        for physical in physical_cartesian_examples() {
            let original_frame = physical.reference_frame();
            let original_mathematical = physical.mathematical_coordinates();
            let mut new_physical = physical.clone();
            
            new_physical.change_reference_frame(target_frame);
            if original_frame != target_frame {
                let changed_mathematical = new_physical.mathematical_coordinates();
                assert!(!original_mathematical.eq_within(changed_mathematical, DISTANCE_ACC));
            }

            new_physical.change_reference_frame(original_frame);
            let mathematical_after_roundtrip = new_physical.mathematical_coordinates();
            assert!(original_mathematical.eq_within(mathematical_after_roundtrip, DISTANCE_ACC));
        }
    }
}
