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
    pub const DISTANCE_ACC: Distance<f64> = Distance { m: ACC };
    pub const ANGLE_ACC: Angle<f64> = Angle { rad: ACC };
}

#[cfg(test)]
pub mod examples {
    use super::*;

    pub fn distance_examples() -> Vec<Distance<f64>> {
        vec![
            Distance::from_m(0.1),
            Distance::from_m(2.0),
            Distance::from_m(-3.0),
        ]
    }

    pub fn positive_distance_examples() -> Vec<Distance<f64>> {
        vec![
            Distance::from_m(0.1),
            Distance::from_m(2.0),
            Distance::from_m(30.),
        ]
    }

    pub fn angle_examples() -> Vec<Angle<f64>> {
        vec![
            Angle::from_degrees(0.0),
            Angle::from_degrees(0.1),
            Angle::from_degrees(45.0),
            Angle::from_degrees(90.0),
            Angle::from_degrees(180.0),
            Angle::from_degrees(270.0),
            Angle::from_degrees(360.0),
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
}

#[cfg(test)]
pub mod benchmarks {
    use super::*;

    pub fn many_angles(num: usize) -> Vec<Angle<f64>> {
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
