use astro_coords::{
    physical_coords::PhysicalCoords,
    reference_frame::{CelestialBody, ReferenceFrame},
    traits::*,
};
use serial_test::serial;
use std::time::{Duration, Instant};

use utils::benchmarks::{many_directions, many_times};

mod utils;

const MAX_DURATION_PER_CALL: Duration = Duration::from_nanos(1000);

#[test]
#[ignore]
#[serial]
fn changing_to_same_frame_is_extremely_fast() {
    const NUM: usize = 1000;
    let dirs = many_directions(NUM);
    let frame = ReferenceFrame::Ecliptic;
    let mut physical_dirs = dirs
        .into_iter()
        .map(|dir| PhysicalCoords::new(dir, frame))
        .collect::<Vec<_>>();

    let start = Instant::now();
    for dir in &mut physical_dirs {
        dir.change_reference_frame(frame);
    }
    let duration = start.elapsed();

    let duration_per_call = duration / NUM as u32;
    println!(
        "{} calls to Direction::change_reference_frame with the same frame took {:?}, or {:?} per call.",
        NUM, duration, duration_per_call
    );
    assert!(duration_per_call < Duration::from_nanos(10));
}

#[test]
#[ignore]
#[serial]
fn changing_from_ecliptic_to_galactic_frame_is_fast() {
    const NUM: usize = 1000;
    let dirs = many_directions(NUM);
    let frame = ReferenceFrame::Ecliptic;
    let mut physical_dirs = dirs
        .into_iter()
        .map(|dir| PhysicalCoords::new(dir, frame))
        .collect::<Vec<_>>();

    let start = Instant::now();
    for dir in &mut physical_dirs {
        dir.change_reference_frame(ReferenceFrame::Galactic);
    }
    let duration = start.elapsed();

    let duration_per_call = duration / NUM as u32;
    println!(
        "{} calls to Direction::change_reference_frame took {:?}, or {:?} per call.",
        NUM, duration, duration_per_call
    );
    assert!(duration_per_call < MAX_DURATION_PER_CALL);
}

#[test]
#[ignore]
#[serial]
fn changing_from_jupiter_to_neptune_cartographic_frame_is_fast() {
    const NUM: usize = 100;
    let dirs = many_directions(NUM);
    let times = many_times(NUM);
    let frame = ReferenceFrame::Cartographic(CelestialBody::Jupiter, times[0]);
    let mut physical_dirs = dirs
        .into_iter()
        .map(|dir| PhysicalCoords::new(dir, frame))
        .collect::<Vec<_>>();
    let total_num = (physical_dirs.len() * times.len()) as u32;

    let start = Instant::now();
    for dir in &mut physical_dirs {
        for time in &times {
            dir.change_reference_frame(ReferenceFrame::Cartographic(CelestialBody::Neptune, *time));
        }
    }
    let duration = start.elapsed();

    let duration_per_call = duration / total_num;
    println!(
        "{} changes from jupiter to nepture cartographic reference frame took {:?}, or {:?} per call.",
        total_num, duration, duration_per_call
    );
    assert!(duration_per_call < MAX_DURATION_PER_CALL);
}
