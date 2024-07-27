use serial_test::serial;
use std::time::{Duration, Instant};

use utils::benchmarks::{many_cartesians, many_directions, many_sphericals};

mod utils;

//TODO

const MAX_DURATION_PER_CALL: Duration = Duration::from_nanos(100);

#[test]
#[ignore]
#[serial]
fn angle_to_for_direction_is_fast() {
    const NUM: usize = 1000;
    let dirs1 = many_directions(NUM);
    let dirs2 = many_directions(NUM);
    let total_calls = NUM * NUM;

    let start = Instant::now();
    let mut dummy: f64 = 0.0;
    for dir1 in &dirs1 {
        for dir2 in &dirs2 {
            dummy += dir1.angle_to(&dir2).rad;
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_call = duration / total_calls as u32;
    println!(
        "{} calls to Direction::angle_to took {:?}, or {:?} per call.",
        total_calls, duration, duration_per_call
    );
    assert!(duration_per_call < MAX_DURATION_PER_CALL);
}

#[test]
#[ignore]
#[serial]
fn angle_to_for_cartesian_is_fast() {
    const NUM: usize = 1000;
    let carts1 = many_cartesians(NUM);
    let carts2 = many_cartesians(NUM);
    let total_calls = NUM * NUM;

    let start = Instant::now();
    let mut dummy: f64 = 0.0;
    for cart1 in &carts1 {
        for cart2 in &carts2 {
            dummy += cart1.angle_to(&cart2).unwrap().rad;
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_call = duration / total_calls as u32;
    println!(
        "{} calls to Cartesian::angle_to took {:?}, or {:?} per call.",
        total_calls, duration, duration_per_call
    );
    assert!(duration_per_call < MAX_DURATION_PER_CALL);
}

#[test]
#[ignore]
#[serial]
fn angle_to_for_spherical_is_fast() {
    const NUM: usize = 1000;
    let sphericals1 = many_sphericals(NUM);
    let sphericals2 = many_sphericals(NUM);
    let total_calls = NUM * NUM;

    let start = Instant::now();
    let mut dummy: f64 = 0.0;
    for sph1 in &sphericals1 {
        for sph2 in &sphericals2 {
            dummy += sph1.angle_to(&sph2).rad;
        }
    }
    let duration = start.elapsed();
    println!("Dummy print: {}", dummy);

    let duration_per_call = duration / total_calls as u32;
    println!(
        "{} calls to Spherical::angle_to took {:?}, or {:?} per call.",
        total_calls, duration, duration_per_call
    );
    assert!(duration_per_call < MAX_DURATION_PER_CALL);
}
