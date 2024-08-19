use astro_coords::traits::*;
use utils::examples::*;

mod utils;

#[test]
fn cartesian_roundtrips() {
    for target_frame in reference_frame_examples() {
        for physical in physical_cartesian_examples() {
            let original_frame = physical.reference_frame();
            let mut new_physical = physical.clone();
            new_physical.change_reference_frame(target_frame);
            if original_frame != target_frame {
                assert_ne!(physical.mathematical_coordinates());
                let length = physical.length();
                let new_cartesian = physical.to_cartesian().unwrap();
                let new_physical = new_cartesian.to_physical_coords(target_frame, length);
                assert!(physical.eq_within(&new_physical, DISTANCE_ACC));
            }
            assert!(cartesian.eq_within(&new_cartesian, DISTANCE_ACC));
        }
    }
}
