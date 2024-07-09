use astro_coords::{direction::Direction, spherical::SphericalCoordinates};
use rand::{rngs::StdRng, Rng, SeedableRng};
use serial_test::serial;
use simple_si_units::geometry::Angle;
use std::{f64::consts::PI, time::Instant};

fn many_directions(num: usize) -> Vec<Direction> {
    let seed = [42; 32];
    let mut rng = StdRng::from_seed(seed);

    let mut directions = Vec::new();
    for _ in 0..num {
        let x = rng.gen_range((-5.)..5.);
        let y = rng.gen_range((-5.)..5.);
        let z = rng.gen_range((-5.)..5.);
        let direction = Direction::new(x, y, z);
        if let Ok(dir) = direction {
            directions.push(dir);
        }
    }
    directions
}

fn many_sphericals(num: usize) -> Vec<SphericalCoordinates> {
    let seed = [42; 32];
    let mut rng = StdRng::from_seed(seed);

    let mut sphericals = Vec::new();
    for _ in 0..num {
        let longitude = Angle {
            rad: rng.gen_range((-PI)..PI),
        };
        let latitude = Angle {
            rad: rng.gen_range((-PI / 2.)..(PI / 2.)),
        };
        sphericals.push(SphericalCoordinates::new(longitude, latitude));
    }
    sphericals
}

#[test]
#[ignore]
#[serial]
fn active_rotation_for_direction_is_fast() {
    const NUM: usize = 1000;
    let vecs = many_directions(NUM);
    let new_zs = many_directions(NUM);
    let total_rotations = vecs.len() * new_zs.len();

    let start = Instant::now();
    let mut dummy = 0.;
    for vec in vecs {
        for new_z in new_zs.iter() {
            let rotated_dir = vec.active_rotation_to_new_z_axis(&new_z);
            dummy += rotated_dir.x();
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_rotation = duration / total_rotations as u32;
    println!(
        "{} calls to Direction::active_rotation_to_new_z_axis took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < std::time::Duration::from_secs(1))
}

#[test]
#[ignore]
#[serial]
fn passive_rotation_for_direction_is_fast() {
    const NUM: usize = 1000;
    let vecs = many_directions(NUM);
    let new_zs = many_directions(NUM);
    let total_rotations = vecs.len() * new_zs.len();

    let start = Instant::now();
    let mut dummy = 0.;
    for vec in vecs {
        for new_z in new_zs.iter() {
            let rotated_dir = vec.passive_rotation_to_new_z_axis(&new_z);
            dummy += rotated_dir.x();
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_rotation = duration / total_rotations as u32;
    println!(
        "{} calls to Direction::passive_rotation_to_new_z_axis took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < std::time::Duration::from_secs(1))
}

#[test]
#[ignore]
#[serial]
fn active_rotation_for_spherical_is_fast() {
    const NUM: usize = 1000;
    let vecs = many_sphericals(NUM);
    let new_zs = many_sphericals(NUM);
    let total_rotations = vecs.len() * new_zs.len();

    let start = Instant::now();
    let mut dummy = 0.;
    for vec in vecs {
        for new_z in new_zs.iter() {
            let rotated_dir = vec.active_rotation_to_new_z_axis(&new_z);
            dummy += rotated_dir.longitude.rad;
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_rotation = duration / total_rotations as u32;
    println!(
        "{} calls to SphericalCoordinates::active_rotation_to_new_z_axis took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < std::time::Duration::from_secs(1))
}

#[test]
#[ignore]
#[serial]
fn passive_rotation_for_spherical_is_fast() {
    const NUM: usize = 1000;
    let vecs = many_sphericals(NUM);
    let new_zs = many_sphericals(NUM);
    let total_rotations = vecs.len() * new_zs.len();

    let start = Instant::now();
    let mut dummy = 0.;
    for vec in vecs {
        for new_z in new_zs.iter() {
            let rotated_dir = vec.passive_rotation_to_new_z_axis(&new_z);
            dummy += rotated_dir.longitude.rad;
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_rotation = duration / total_rotations as u32;
    println!(
        "{} calls to SphericalCoordinates::passive_rotation_to_new_z_axis took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < std::time::Duration::from_secs(1))
}
