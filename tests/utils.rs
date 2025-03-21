#![allow(dead_code)]

use std::f64::consts::PI;

use rand::{Rng, SeedableRng, rngs::StdRng};

use astro_coords::{
    cartesian::Cartesian, direction::Direction, earth_equatorial::EarthEquatorial,
    ecliptic::Ecliptic, equatorial::Equatorial, spherical::Spherical,
};

#[cfg(test)]
pub mod constants {
    use uom::si::{
        angle::radian,
        f64::{Angle, Length},
        length::meter,
    };

    pub const ACC: f64 = 1e-5;
    pub fn distance_acc() -> Length {
        Length::new::<meter>(ACC)
    }
    pub fn angle_acc() -> Angle {
        Angle::new::<radian>(ACC)
    }
}

#[cfg(test)]
pub mod examples {
    use astro_coords::{
        physical_coords::PhysicalCoords,
        reference_frame::{CelestialBody, ReferenceFrame, RotationalElements},
    };
    use uom::si::{
        angle::degree,
        f64::{Angle, Length, Time},
        length::meter,
        time::{day, year},
    };

    use super::*;

    pub fn distance_examples() -> Vec<Length> {
        vec![
            Length::new::<meter>(0.1),
            Length::new::<meter>(2.0),
            Length::new::<meter>(-3.0),
        ]
    }

    pub fn positive_distance_examples() -> Vec<Length> {
        vec![
            Length::new::<meter>(0.1),
            Length::new::<meter>(2.0),
            Length::new::<meter>(30.),
        ]
    }

    pub fn angle_examples() -> Vec<Angle> {
        vec![
            Angle::new::<degree>(0.0),
            Angle::new::<degree>(0.1),
            Angle::new::<degree>(45.0),
            Angle::new::<degree>(90.0),
            Angle::new::<degree>(180.0),
            Angle::new::<degree>(270.0),
            Angle::new::<degree>(360.0),
        ]
    }

    pub fn long_time_examples() -> Vec<Time> {
        vec![
            Time::new::<year>(0.0),
            Time::new::<year>(1.0),
            Time::new::<year>(100.0),
            Time::new::<year>(1000.0),
            Time::new::<year>(10000.0),
        ]
    }

    pub fn spherical_examples() -> Vec<Spherical> {
        let mut examples = vec![];
        for lon in angle_examples() {
            for lat in angle_examples() {
                examples.push(Spherical::new(lon, lat));
            }
        }
        examples
    }

    pub fn cartesian_examples() -> Vec<Cartesian> {
        let mut examples = vec![];
        for x in distance_examples() {
            for y in distance_examples() {
                for z in distance_examples() {
                    examples.push(Cartesian::new(x, y, z));
                }
            }
        }
        examples
    }

    pub fn direction_examples() -> Vec<Direction> {
        let mut examples = vec![];
        for x in distance_examples() {
            for y in distance_examples() {
                for z in distance_examples() {
                    let dir = Direction::new(x.value, y.value, z.value).unwrap();
                    examples.push(dir);
                }
            }
        }
        examples
    }

    pub fn earth_equatorial_examples() -> Vec<EarthEquatorial> {
        let mut examples = vec![];
        for ra in angle_examples() {
            for dec in angle_examples() {
                examples.push(EarthEquatorial::new(ra, dec));
            }
        }
        examples
    }

    pub fn ecliptic_examples() -> Vec<Ecliptic> {
        let mut examples = vec![];
        for spherical in spherical_examples() {
            examples.push(Ecliptic::new(spherical));
        }
        examples
    }

    pub fn equatorial_examples() -> Vec<Equatorial> {
        let mut examples = vec![];
        for spherical in spherical_examples() {
            for axis in direction_examples() {
                examples.push(Equatorial::new(spherical, axis));
            }
        }
        examples
    }

    pub fn celestial_body_examples() -> Vec<CelestialBody> {
        let mut vec = vec![
            CelestialBody::Sun,
            CelestialBody::Mercury,
            CelestialBody::Venus,
            CelestialBody::Earth,
            CelestialBody::Mars,
            CelestialBody::Jupiter,
            CelestialBody::Saturn,
            CelestialBody::Uranus,
            CelestialBody::Neptune,
        ];
        let sphericals = spherical_examples();
        let angles = angle_examples();
        let angular_velocities: Vec<_> = angle_examples()
            .into_iter()
            .map(|a| a / Time::new::<day>(1.0))
            .collect();
        for i in 0..sphericals.len() {
            let j = i * 12345 % angles.len();
            let k = i * 54321 % angular_velocities.len();
            let rot = RotationalElements {
                z_axis: sphericals[i],
                prime_meridian_offset_offset: angles[j],
                prime_meridian_offset_rate: angular_velocities[k].into(),
            };
            vec.push(CelestialBody::Custom(rot));
        }
        vec
    }

    pub fn reference_frame_examples() -> Vec<ReferenceFrame> {
        let mut vec = vec![
            ReferenceFrame::Equatorial,
            ReferenceFrame::Ecliptic,
            ReferenceFrame::Galactic,
        ];
        for body in celestial_body_examples() {
            for time in long_time_examples() {
                vec.push(ReferenceFrame::Cartographic(body, time));
            }
        }
        vec
    }

    pub fn physical_cartesian_examples() -> Vec<PhysicalCoords<Cartesian>> {
        let mut examples = vec![];
        for frame in reference_frame_examples() {
            for cartesian in cartesian_examples() {
                examples.push(PhysicalCoords::new(cartesian, frame));
            }
        }
        examples
    }

    pub fn physical_direction_examples() -> Vec<PhysicalCoords<Direction>> {
        let mut examples = vec![];
        for frame in reference_frame_examples() {
            for direction in direction_examples() {
                examples.push(PhysicalCoords::new(direction, frame));
            }
        }
        examples
    }

    pub fn physical_spherical_examples() -> Vec<PhysicalCoords<Spherical>> {
        let mut examples = vec![];
        for frame in reference_frame_examples() {
            for spherical in spherical_examples() {
                examples.push(PhysicalCoords::new(spherical, frame));
            }
        }
        examples
    }
}

#[cfg(test)]
pub mod benchmarks {

    use uom::si::{
        angle::radian,
        f64::{Angle, Length, Time},
        length::meter,
        time::year,
    };

    use super::*;

    pub fn many_angles(num: usize) -> Vec<Angle> {
        let seed = [42; 32];
        let mut rng = StdRng::from_seed(seed);

        let mut angles = Vec::new();
        for _ in 0..num {
            let rad = rng.random_range((-PI)..PI);
            angles.push(Angle::new::<radian>(rad));
        }
        angles
    }

    pub fn many_directions(num: usize) -> Vec<Direction> {
        let seed = [42; 32];
        let mut rng = StdRng::from_seed(seed);

        let mut directions = Vec::new();
        for _ in 0..num {
            let x = rng.random_range((-5.)..5.);
            let y = rng.random_range((-5.)..5.);
            let z = rng.random_range((-5.)..5.);
            let direction = Direction::new(x, y, z);
            if let Ok(dir) = direction {
                directions.push(dir);
            }
        }
        directions
    }

    pub fn many_times(num: usize) -> Vec<Time> {
        let seed = [42; 32];
        let mut rng = StdRng::from_seed(seed);

        let mut times = Vec::new();
        for _ in 0..num {
            let yr = rng.random_range(0.0..1000.0);
            times.push(Time::new::<year>(yr));
        }
        times
    }

    pub fn many_cartesians(num: usize) -> Vec<Cartesian> {
        let seed = [42; 32];
        let mut rng = StdRng::from_seed(seed);

        let mut cartesians = Vec::new();
        for _ in 0..num {
            let x = rng.random_range((-5.)..5.);
            let y = rng.random_range((-5.)..5.);
            let z = rng.random_range((-5.)..5.);
            cartesians.push(Cartesian::new(
                Length::new::<meter>(x),
                Length::new::<meter>(y),
                Length::new::<meter>(z),
            ));
        }
        cartesians
    }

    pub fn many_sphericals(num: usize) -> Vec<Spherical> {
        let seed = [42; 32];
        let mut rng = StdRng::from_seed(seed);

        let mut sphericals = Vec::new();
        for _ in 0..num {
            let longitude = Angle::new::<radian>(rng.random_range((-PI)..PI));
            let latitude = Angle::new::<radian>(rng.random_range((-PI / 2.)..(PI / 2.)));
            sphericals.push(Spherical::new(longitude, latitude));
        }
        sphericals
    }
}
