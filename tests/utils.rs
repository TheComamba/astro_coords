use std::f64::consts::PI;

use rand::{rngs::StdRng, Rng, SeedableRng};
use simple_si_units::{base::Distance, geometry::Angle};

use astro_coords::{
    cartesian::Cartesian, direction::Direction, earth_equatorial::EarthEquatorial,
    ecliptic::Ecliptic, equatorial::Equatorial, spherical::Spherical,
};

#[cfg(test)]
pub mod constants {
    use super::*;

    pub const ACC: f64 = 1e-5;
    pub const DISTANCE_ACC: Length = Distance { m: ACC };
    pub const ANGLE_ACC: Angle = Angle { rad: ACC };
}

#[cfg(test)]
pub mod examples {
    use astro_coords::{
        physical_coords::PhysicalCoords,
        reference_frame::{CelestialBody, ReferenceFrame, RotationalElements},
    };
    use simple_si_units::base::Time;

    use super::*;

    pub fn distance_examples() -> Vec<Length> {
        vec![
            Distance::from_m(0.1),
            Distance::from_m(2.0),
            Distance::from_m(-3.0),
        ]
    }

    pub fn positive_distance_examples() -> Vec<Length> {
        vec![
            Distance::from_m(0.1),
            Distance::from_m(2.0),
            Distance::from_m(30.),
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
            Time::from_yr(0.0),
            Time::from_yr(1.0),
            Time::from_yr(100.0),
            Time::from_yr(1000.0),
            Time::from_yr(10000.0),
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
                    let dir = Direction::new(x.m, y.m, z.m).unwrap();
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
            .map(|a| a / Time::from_days(1.0))
            .collect();
        for i in 0..sphericals.len() {
            let j = i * 12345 % angles.len();
            let k = i * 54321 % angular_velocities.len();
            let rot = RotationalElements {
                z_axis: sphericals[i],
                prime_meridian_offset_offset: angles[j],
                prime_meridian_offset_rate: angular_velocities[k],
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
    use simple_si_units::base::Time;

    use super::*;

    pub fn many_angles(num: usize) -> Vec<Angle> {
        let seed = [42; 32];
        let mut rng = StdRng::from_seed(seed);

        let mut angles = Vec::new();
        for _ in 0..num {
            let rad = rng.gen_range((-PI)..PI);
            angles.push(Angle { rad });
        }
        angles
    }

    pub fn many_directions(num: usize) -> Vec<Direction> {
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

    pub fn many_times(num: usize) -> Vec<Time> {
        let seed = [42; 32];
        let mut rng = StdRng::from_seed(seed);

        let mut times = Vec::new();
        for _ in 0..num {
            let yr = rng.gen_range(0.0..1000.0);
            times.push(Time::from_yr(yr));
        }
        times
    }

    pub fn many_cartesians(num: usize) -> Vec<Cartesian> {
        let seed = [42; 32];
        let mut rng = StdRng::from_seed(seed);

        let mut cartesians = Vec::new();
        for _ in 0..num {
            let x = rng.gen_range((-5.)..5.);
            let y = rng.gen_range((-5.)..5.);
            let z = rng.gen_range((-5.)..5.);
            cartesians.push(Cartesian::new(
                Distance::from_m(x),
                Distance::from_m(y),
                Distance::from_m(z),
            ));
        }
        cartesians
    }

    pub fn many_sphericals(num: usize) -> Vec<Spherical> {
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
            sphericals.push(Spherical::new(longitude, latitude));
        }
        sphericals
    }
}
