use simple_si_units::base::Distance;
use utils::{constants::*, examples::*};

mod utils;

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
        let ecliptic = cartesian.to_ecliptic().unwrap();
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
        let spherical = cartesian.to_spherical().unwrap();
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
        let new_ecliptic = cartesian.to_ecliptic().unwrap();
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
        let axis = equatorial.rotation_axis.clone();
        let length = Distance::from_m(10.);
        let cartesian = equatorial.to_cartesian(length);
        let new_equatorial = cartesian.to_equatorial(axis).unwrap();
        assert!(equatorial.eq_within(&new_equatorial, ANGLE_ACC));
    }
}

#[test]
fn equatorial_to_direction_roundtrip() {
    for equatorial in equatorial_examples() {
        let axis = equatorial.rotation_axis.clone();
        let direction = equatorial.to_direction();
        let new_equatorial = direction.to_equatorial(axis);
        assert!(equatorial.eq_within(&new_equatorial, ANGLE_ACC));
    }
}

#[test]
fn equatorial_to_earth_equatorial_roundtrip() {
    for equatorial in equatorial_examples() {
        let axis = equatorial.rotation_axis.clone();
        let earth_equatorial = equatorial.to_earth_equatorial();
        let new_equatorial = earth_equatorial.to_equatorial(axis);
        assert!(equatorial.eq_within(&new_equatorial, ANGLE_ACC));
    }
}

#[test]
fn equatorial_to_ecliptic_roundtrip() {
    for equatorial in equatorial_examples() {
        let axis = equatorial.rotation_axis.clone();
        let ecliptic = equatorial.to_ecliptic();
        let new_equatorial = ecliptic.to_equatorial(axis);
        assert!(equatorial.eq_within(&new_equatorial, ANGLE_ACC));
    }
}

#[test]
fn equatorial_to_spherical_roundtrip() {
    for equatorial in equatorial_examples() {
        let axis = equatorial.rotation_axis.clone();
        let spherical = equatorial.to_spherical();
        let new_equatorial = spherical.to_equatorial(axis);
        assert!(equatorial.eq_within(&new_equatorial, ANGLE_ACC));
    }
}

#[test]
fn spherical_to_cartesian_roundtrip() {
    for spherical in spherical_examples() {
        let length = Distance::from_m(10.);
        let cartesian = spherical.to_cartesian(length);
        let new_spherical = cartesian.to_spherical().unwrap();
        assert!(spherical.eq_within(&new_spherical, ANGLE_ACC));
    }
}

#[test]
fn spherical_to_direction_roundtrip() {
    for spherical in spherical_examples() {
        let direction = spherical.to_direction();
        let new_spherical = direction.to_spherical();
        assert!(spherical.eq_within(&new_spherical, ANGLE_ACC));
    }
}

#[test]
fn spherical_to_earth_equatorial_roundtrip() {
    for spherical in spherical_examples() {
        let earth_equatorial = spherical.to_earth_equatorial();
        let new_spherical = earth_equatorial.to_spherical();
        assert!(spherical.eq_within(&new_spherical, ANGLE_ACC));
    }
}

#[test]
fn spherical_to_ecliptic_roundtrip() {
    for spherical in spherical_examples() {
        let ecliptic = spherical.to_ecliptic();
        let new_spherical = ecliptic.to_spherical();
        assert!(spherical.eq_within(&new_spherical, ANGLE_ACC));
    }
}

#[test]
fn spherical_to_equatorial_roundtrip() {
    for spherical in spherical_examples() {
        for axis in direction_examples() {
            let equatorial = spherical.to_equatorial(axis);
            let new_spherical = equatorial.to_spherical();
            assert!(spherical.eq_within(&new_spherical, ANGLE_ACC));
        }
    }
}
