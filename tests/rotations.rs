use astro_coords::direction::Direction;
use simple_si_units::geometry::Angle;

const ACC: Angle<f64> = Angle { rad: 1e-5 };

fn example_directions() -> Vec<Direction> {
    let ordinates = vec![0.0, 0.5, -0.1, 12.0];
    let mut directions = Vec::new();
    for x in ordinates.iter() {
        for y in ordinates.iter() {
            for z in ordinates.iter() {
                let direction = Direction::new(*x, *y, *z);
                if let Ok(dir) = direction {
                    directions.push(dir);
                }
            }
        }
    }
    directions
}

#[test]
fn active_rotation_is_the_same_for_direction_and_spherical() {
    for new_z in example_directions() {
        for vec in example_directions() {
            let rotated_dir = vec.active_rotation_to_new_z_axis(&new_z);

            let spherical_new_z = new_z.to_spherical();
            let spherical_dir = vec.to_spherical();
            let rotated_spherical = spherical_dir.active_rotation_to_new_z_axis(&spherical_new_z);

            let expected = rotated_dir.to_spherical();
            assert!(
                rotated_spherical.eq_within(&expected, ACC),
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
    for new_z in example_directions() {
        for vec in example_directions() {
            let rotated_dir = vec.passive_rotation_to_new_z_axis(&new_z);

            let spherical_new_z = new_z.to_spherical();
            let spherical_dir = vec.to_spherical();
            let rotated_spherical = spherical_dir.passive_rotation_to_new_z_axis(&spherical_new_z);

            let expected = rotated_dir.to_spherical();
            assert!(rotated_spherical.eq_within(&expected, ACC));
        }
    }
}
