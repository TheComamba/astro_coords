use astro_coords::traits::*;
use serial_test::serial;
use std::time::{Duration, Instant};

use utils::benchmarks::*;

mod utils;

const MAX_DURATION_PER_SIMPLE_ROTATION_FOR_DIRECTION: Duration = Duration::from_nanos(100);
const MAX_DURATION_PER_SIMPLE_ROTATION_FOR_SPHERICAL: Duration = Duration::from_nanos(1000);
const MAX_DURATION_PER_ELABORATE_ROTATION: Duration = Duration::from_nanos(1000);

#[test]
#[ignore]
#[serial]
fn rotation_for_direction_is_fast() {
    const NUM: usize = 100;
    let angles = many_angles(NUM);
    let vecs = many_directions(NUM);
    let axes = many_directions(NUM);
    let total_rotations = angles.len() * vecs.len() * axes.len();

    let start = Instant::now();
    let mut dummy = 0.;
    for angle in angles {
        for vec in vecs.iter() {
            for axis in axes.iter() {
                let rotated_dir = vec.rotated(angle, &axis);
                dummy += rotated_dir.x();
            }
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_rotation = duration / total_rotations as u32;
    println!(
        "{} calls to Direction::rotated took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < MAX_DURATION_PER_SIMPLE_ROTATION_FOR_DIRECTION)
}

#[test]
#[ignore]
#[serial]
fn x_rotation_for_direction_is_fast() {
    const NUM: usize = 1000;
    let angles = many_angles(NUM);
    let vecs = many_directions(NUM);
    let total_rotations = angles.len() * vecs.len();

    let start = Instant::now();
    let mut dummy = 0.;
    for angle in angles {
        for vec in vecs.iter() {
            let rotated_dir = vec.rotated_x(angle);
            dummy += rotated_dir.y();
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_rotation = duration / total_rotations as u32;
    println!(
        "{} calls to Direction::rotated_x took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < MAX_DURATION_PER_SIMPLE_ROTATION_FOR_DIRECTION)
}

#[test]
#[ignore]
#[serial]
fn y_rotation_for_direction_is_fast() {
    const NUM: usize = 1000;
    let angles = many_angles(NUM);
    let vecs = many_directions(NUM);
    let total_rotations = angles.len() * vecs.len();

    let start = Instant::now();
    let mut dummy = 0.;
    for angle in angles {
        for vec in vecs.iter() {
            let rotated_dir = vec.rotated_y(angle);
            dummy += rotated_dir.x();
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_rotation = duration / total_rotations as u32;
    println!(
        "{} calls to Direction::rotated_y took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < MAX_DURATION_PER_SIMPLE_ROTATION_FOR_DIRECTION)
}

#[test]
#[ignore]
#[serial]
fn z_rotation_for_direction_is_fast() {
    const NUM: usize = 1000;
    let angles = many_angles(NUM);
    let vecs = many_directions(NUM);
    let total_rotations = angles.len() * vecs.len();

    let start = Instant::now();
    let mut dummy = 0.;
    for angle in angles {
        for vec in vecs.iter() {
            let rotated_dir = vec.rotated_z(angle);
            dummy += rotated_dir.x();
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_rotation = duration / total_rotations as u32;
    println!(
        "{} calls to Direction::rotated_z took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < MAX_DURATION_PER_SIMPLE_ROTATION_FOR_DIRECTION)
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
    assert!(duration_per_rotation < MAX_DURATION_PER_ELABORATE_ROTATION)
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
    assert!(duration_per_rotation < MAX_DURATION_PER_ELABORATE_ROTATION)
}

#[test]
#[ignore]
#[serial]
fn rotation_for_sphericals_is_fast() {
    const NUM: usize = 100;
    let angles = many_angles(NUM);
    let vecs = many_sphericals(NUM);
    let axes = many_directions(NUM);
    let total_rotations = angles.len() * vecs.len() * axes.len();

    let start = Instant::now();
    let mut dummy = 0.;
    for angle in angles {
        for vec in vecs.iter() {
            for axis in axes.iter() {
                let rotated_dir = vec.rotated(angle, &axis);
                dummy += rotated_dir.longitude.rad;
            }
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_rotation = duration / total_rotations as u32;
    println!(
        "{} calls to Spherical::rotated took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < MAX_DURATION_PER_SIMPLE_ROTATION_FOR_SPHERICAL)
}

#[test]
#[ignore]
#[serial]
fn x_rotation_for_sphericals_is_fast() {
    const NUM: usize = 1000;
    let angles = many_angles(NUM);
    let vecs = many_sphericals(NUM);
    let total_rotations = angles.len() * vecs.len();

    let start = Instant::now();
    let mut dummy = 0.;
    for angle in angles {
        for vec in vecs.iter() {
            let rotated_dir = vec.rotated_x(angle);
            dummy += rotated_dir.longitude.rad;
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_rotation = duration / total_rotations as u32;
    println!(
        "{} calls to Spherical::rotated_x took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < MAX_DURATION_PER_SIMPLE_ROTATION_FOR_SPHERICAL)
}

#[test]
#[ignore]
#[serial]
fn y_rotation_for_sphericals_is_fast() {
    const NUM: usize = 1000;
    let angles = many_angles(NUM);
    let vecs = many_sphericals(NUM);
    let total_rotations = angles.len() * vecs.len();

    let start = Instant::now();
    let mut dummy = 0.;
    for angle in angles {
        for vec in vecs.iter() {
            let rotated_dir = vec.rotated_y(angle);
            dummy += rotated_dir.longitude.rad;
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_rotation = duration / total_rotations as u32;
    println!(
        "{} calls to Spherical::rotated_y took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < MAX_DURATION_PER_SIMPLE_ROTATION_FOR_SPHERICAL)
}

#[test]
#[ignore]
#[serial]
fn z_rotation_for_sphericals_is_fast() {
    const NUM: usize = 1000;
    let angles = many_angles(NUM);
    let vecs = many_sphericals(NUM);
    let total_rotations = angles.len() * vecs.len();

    let start = Instant::now();
    let mut dummy = 0.;
    for angle in angles {
        for vec in vecs.iter() {
            let rotated_dir = vec.rotated_z(angle);
            dummy += rotated_dir.longitude.rad;
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_rotation = duration / total_rotations as u32;
    println!(
        "{} calls to Spherical::rotated_z took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < MAX_DURATION_PER_SIMPLE_ROTATION_FOR_SPHERICAL)
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
        "{} calls to Spherical::active_rotation_to_new_z_axis took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < MAX_DURATION_PER_ELABORATE_ROTATION)
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
        "{} calls to Spherical::passive_rotation_to_new_z_axis took {:?}, or {:?} per call.",
        total_rotations, duration, duration_per_rotation
    );
    assert!(duration_per_rotation < MAX_DURATION_PER_ELABORATE_ROTATION)
}
