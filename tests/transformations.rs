use astro_coordinates::{
    cartesian::CartesianCoordinates,
    declination::{Declination, Sgn},
    direction::Direction,
    earth_equatorial::EarthEquatorialCoordinates,
    ecliptic::EclipticCoordinates,
    equatorial::EquatorialCoordinates,
    right_ascension::RightAscension,
    spherical::SphericalCoordinates,
};
use simple_si_units::{base::Distance, geometry::Angle};

const ACC: f64 = 1e-5;
const DISTANCE_ACC: Distance<f64> = Distance { m: ACC };
const ANGLE_ACC: Angle<f64> = Angle { rad: ACC };

fn distance_examples() -> Vec<Distance<f64>> {
    vec![
        Distance::from_m(0.1),
        Distance::from_m(2.0),
        Distance::from_m(-3.0),
    ]
}

fn angle_examples() -> Vec<Angle<f64>> {
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

fn spherical_examples() -> Vec<SphericalCoordinates> {
    let mut examples = vec![];
    for lon in angle_examples() {
        for lat in angle_examples() {
            examples.push(SphericalCoordinates::new(lon, lat));
        }
    }
    examples
}

fn declination_examples() -> Vec<Declination> {
    vec![
        Declination::new(Sgn::Pos, 0, 0, 0),
        Declination::new(Sgn::Neg, 0, 0, 0),
        Declination::new(Sgn::Pos, 90, 0, 0),
        Declination::new(Sgn::Neg, 90, 0, 0),
        Declination::new(Sgn::Pos, 45, 0, 0),
        Declination::new(Sgn::Neg, 45, 0, 0),
        Declination::new(Sgn::Pos, 0, 30, 0),
        Declination::new(Sgn::Neg, 0, 30, 0),
        Declination::new(Sgn::Pos, 0, 0, 30),
        Declination::new(Sgn::Neg, 0, 0, 30),
        Declination::new(Sgn::Pos, 45, 30, 30),
        Declination::new(Sgn::Neg, 45, 30, 30),
        Declination::new(Sgn::Pos, 90, 30, 30),
        Declination::new(Sgn::Neg, 90, 30, 30),
    ]
}

fn ra_examples() -> Vec<RightAscension> {
    vec![
        RightAscension::new(0, 0, 0),
        RightAscension::new(6, 0, 0),
        RightAscension::new(12, 0, 0),
        RightAscension::new(18, 0, 0),
        RightAscension::new(23, 59, 59),
        RightAscension::new(0, 59, 59),
        RightAscension::new(6, 59, 59),
        RightAscension::new(12, 59, 59),
        RightAscension::new(18, 59, 59),
    ]
}

fn cartesian_examples() -> Vec<CartesianCoordinates> {
    let mut examples = vec![];
    for x in distance_examples() {
        for y in distance_examples() {
            for z in distance_examples() {
                examples.push(CartesianCoordinates::new(x, y, z));
            }
        }
    }
    examples
}

fn direction_examples() -> Vec<Direction> {
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

fn earth_equatorial_examples() -> Vec<EarthEquatorialCoordinates> {
    let mut examples = vec![];
    for ra in angle_examples() {
        for dec in angle_examples() {
            examples.push(EarthEquatorialCoordinates::new(ra, dec));
        }
    }
    examples
}

fn ecliptic_examples() -> Vec<EclipticCoordinates> {
    let mut examples = vec![];
    for spherical in spherical_examples() {
        examples.push(EclipticCoordinates::new(spherical));
    }
    examples
}

fn equatorial_examples() -> Vec<EquatorialCoordinates> {
    let mut examples = vec![];
    for spherical in spherical_examples() {
        for axis in direction_examples() {
            examples.push(EquatorialCoordinates::new(spherical, axis));
        }
    }
    examples
}

#[test]
fn cartesian_to_direction_roundtrip() {
    for cartesian in cartesian_examples() {
        let length = cartesian.length();
        let direction = cartesian.to_direction().unwrap();
        let new_cartesian = direction.to_cartesian(length);
        assert!(cartesian.eq_within(&new_cartesian, DISTANCE_ACC));
    }
}

#[test]
fn cartesian_to_earth_equatorial_roundtrip() {
    for cartesian in cartesian_examples() {
        let length = cartesian.length();
        let earth_equatorial = cartesian.to_earth_equatorial().unwrap();
        let new_cartesian = earth_equatorial.to_cartesian(length);
        assert!(cartesian.eq_within(&new_cartesian, DISTANCE_ACC));
    }
}

#[test]
fn cartesian_to_ecliptic_roundtrip() {
    for cartesian in cartesian_examples() {
        let length = cartesian.length();
        let ecliptic = cartesian.to_ecliptic();
        let new_cartesian = ecliptic.to_cartesian(length);
        assert!(cartesian.eq_within(&new_cartesian, DISTANCE_ACC));
    }
}

#[test]
fn cartesian_to_equatorial_roundtrip() {
    for cartesian in cartesian_examples() {
        for axis in direction_examples() {
            let length = cartesian.length();
            let equatorial = cartesian.to_equatorial(axis).unwrap();
            let new_cartesian = equatorial.to_cartesian(length);
            assert!(cartesian.eq_within(&new_cartesian, DISTANCE_ACC));
        }
    }
}

#[test]
fn cartesian_to_spherical_roundtrip() {
    for cartesian in cartesian_examples() {
        let length = cartesian.length();
        let spherical = cartesian.to_spherical();
        let new_cartesian = spherical.to_cartesian(length);
        assert!(cartesian.eq_within(&new_cartesian, DISTANCE_ACC));
    }
}

#[test]
fn direction_to_cartesian_roundtrip() {
    for direction in direction_examples() {
        let length = Distance::from_m(10.);
        let cartesian = direction.to_cartesian(length);
        let new_direction = cartesian.to_direction().unwrap();
        assert!(direction.eq_within(&new_direction, ACC));
    }
}

#[test]
fn direction_to_earth_equatorial_roundtrip() {
    for direction in direction_examples() {
        let earth_equatorial = direction.to_earth_equatorial();
        let new_direction = earth_equatorial.to_direction();
        assert!(direction.eq_within(&new_direction, ACC));
    }
}

#[test]
fn direction_to_ecliptic_roundtrip() {
    for direction in direction_examples() {
        let ecliptic = direction.to_ecliptic();
        let new_direction = ecliptic.to_direction();
        assert!(direction.eq_within(&new_direction, ACC));
    }
}

#[test]
fn direction_to_equatorial_roundtrip() {
    for direction in direction_examples() {
        for axis in direction_examples() {
            let equatorial = direction.to_equatorial(axis);
            let new_direction = equatorial.to_direction();
            assert!(direction.eq_within(&new_direction, ACC));
        }
    }
}

#[test]
fn direction_to_spherical_roundtrip() {
    for direction in direction_examples() {
        let spherical = direction.to_spherical();
        let new_direction = spherical.to_direction();
        assert!(direction.eq_within(&new_direction, ACC));
    }
}

#[test]
fn earth_equatorial_to_cartesian_roundtrip() {
    for earth_equatorial in earth_equatorial_examples() {
        let length = Distance::from_m(10.);
        let cartesian = earth_equatorial.to_cartesian(length);
        let new_earth_equatorial = cartesian.to_earth_equatorial().unwrap();
        assert!(earth_equatorial.eq_within(&new_earth_equatorial, ANGLE_ACC));
    }
}

#[test]
fn earth_equatorial_to_direction_roundtrip() {
    for earth_equatorial in earth_equatorial_examples() {
        let direction = earth_equatorial.to_direction();
        let new_earth_equatorial = direction.to_earth_equatorial();
        assert!(earth_equatorial.eq_within(&new_earth_equatorial, ANGLE_ACC));
    }
}

#[test]
fn earth_equatorial_to_ecliptic_roundtrip() {
    for earth_equatorial in earth_equatorial_examples() {
        let ecliptic = earth_equatorial.to_ecliptic();
        let new_earth_equatorial = ecliptic.to_earth_equatorial();
        assert!(earth_equatorial.eq_within(&new_earth_equatorial, ANGLE_ACC));
    }
}

#[test]
fn earth_equatorial_to_equatorial_roundtrip() {
    for earth_equatorial in earth_equatorial_examples() {
        for axis in direction_examples() {
            let equatorial = earth_equatorial.to_equatorial(axis);
            let new_earth_equatorial = equatorial.to_earth_equatorial();
            assert!(earth_equatorial.eq_within(&new_earth_equatorial, ANGLE_ACC));
        }
    }
}

#[test]
fn earth_equatorial_to_spherical_roundtrip() {
    for earth_equatorial in earth_equatorial_examples() {
        let spherical = earth_equatorial.to_spherical();
        let new_earth_equatorial = spherical.to_earth_equatorial();
        assert!(earth_equatorial.eq_within(&new_earth_equatorial, ANGLE_ACC));
    }
}

#[test]
fn ecliptic_to_cartesian_roundtrip() {
    for ecliptic in ecliptic_examples() {
        let length = Distance::from_m(10.);
        let cartesian = ecliptic.to_cartesian(length);
        let new_ecliptic = cartesian.to_ecliptic();
        assert!(ecliptic.eq_within(&new_ecliptic, ANGLE_ACC));
    }
}

#[test]
fn ecliptic_to_direction_roundtrip() {
    for ecliptic in ecliptic_examples() {
        let direction = ecliptic.to_direction();
        let new_ecliptic = direction.to_ecliptic();
        assert!(ecliptic.eq_within(&new_ecliptic, ANGLE_ACC));
    }
}

#[test]
fn ecliptic_to_earth_equatorial_roundtrip() {
    for ecliptic in ecliptic_examples() {
        let earth_equatorial = ecliptic.to_earth_equatorial();
        let new_ecliptic = earth_equatorial.to_ecliptic();
        assert!(ecliptic.eq_within(&new_ecliptic, ANGLE_ACC));
    }
}

#[test]
fn ecliptic_to_equatorial_roundtrip() {
    for ecliptic in ecliptic_examples() {
        for axis in direction_examples() {
            let equatorial = ecliptic.to_equatorial(axis);
            let new_ecliptic = equatorial.to_ecliptic();
            assert!(ecliptic.eq_within(&new_ecliptic, ANGLE_ACC));
        }
    }
}

#[test]
fn ecliptic_to_spherical_roundtrip() {
    for ecliptic in ecliptic_examples() {
        let spherical = ecliptic.to_spherical();
        let new_ecliptic = spherical.to_ecliptic();
        assert!(ecliptic.eq_within(&new_ecliptic, ANGLE_ACC));
    }
}

#[test]
fn equatorial_to_cartesian_roundtrip() {
    for equatorial in equatorial_examples() {
        let axis = equatorial.get_rotation_axis().clone();
        let length = Distance::from_m(10.);
        let cartesian = equatorial.to_cartesian(length);
        let new_equatorial = cartesian.to_equatorial(axis).unwrap();
        assert!(equatorial.eq_within(&new_equatorial, ANGLE_ACC));
    }
}

#[test]
fn equatorial_to_direction_roundtrip() {
    for equatorial in equatorial_examples() {
        let axis = equatorial.get_rotation_axis().clone();
        let direction = equatorial.to_direction();
        let new_equatorial = direction.to_equatorial(axis);
        assert!(equatorial.eq_within(&new_equatorial, ANGLE_ACC));
    }
}

#[test]
fn equatorial_to_earth_equatorial_roundtrip() {
    for equatorial in equatorial_examples() {
        let axis = equatorial.get_rotation_axis().clone();
        let earth_equatorial = equatorial.to_earth_equatorial();
        let new_equatorial = earth_equatorial.to_equatorial(axis);
        assert!(equatorial.eq_within(&new_equatorial, ANGLE_ACC));
    }
}

#[test]
fn equatorial_to_ecliptic_roundtrip() {
    for equatorial in equatorial_examples() {
        let axis = equatorial.get_rotation_axis().clone();
        let ecliptic = equatorial.to_ecliptic();
        let new_equatorial = ecliptic.to_equatorial(axis);
        assert!(equatorial.eq_within(&new_equatorial, ANGLE_ACC));
    }
}

#[test]
fn equatorial_to_spherical_roundtrip() {
    for equatorial in equatorial_examples() {
        let axis = equatorial.get_rotation_axis().clone();
        let spherical = equatorial.to_spherical();
        let new_equatorial = spherical.to_equatorial(axis);
        assert!(equatorial.eq_within(&new_equatorial, ANGLE_ACC));
    }
}
