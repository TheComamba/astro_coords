use std::f64::consts::PI;

use uom::si::{
    angle::{degree, radian},
    f64::Angle,
};

pub(crate) fn ANGLE_ZERO() -> Angle {
    Angle::new::<radian>(0.)
}
pub(crate) fn FULL_CIRC() -> Angle {
    Angle::new::<radian>(2. * PI)
}
pub(crate) fn QUARTER_CIRC() -> Angle {
    Angle::new::<radian>(2. * PI / 4.)
}
pub(crate) fn HALF_CIRC() -> Angle {
    Angle::new::<radian>(2. * PI / 2.)
}
pub(crate) fn EARTH_AXIS_TILT() -> Angle {
    Angle::new::<degree>(23.439_281)
}
pub(crate) fn GALACTIC_LONGITUDE_OF_NORTH_CELESTIAL_POLE() -> Angle {
    Angle::new::<radian>(122.93314)
}
pub(crate) fn RA_OF_GALACTIC_NORTH() -> Angle {
    Angle::new::<radian>(3.3658674624710652)
}
pub(crate) fn DEC_OF_GALACTIC_NORTH() -> Angle {
    Angle::new::<radian>(27.13)
}

pub(crate) fn angle_eq_within(actual: Angle, expected: Angle, accuracy: Angle) -> bool {
    let diff = normalized_angle(actual - expected);
    diff.get::<radian>().abs() < accuracy.get::<radian>()
}

/// Normalize the angle to a range of −π to +π radians, -180° to 180°.
pub fn normalized_angle(mut angle: Angle) -> Angle {
    angle %= FULL_CIRC();
    if angle > HALF_CIRC() {
        angle -= FULL_CIRC();
    } else if angle < -HALF_CIRC() {
        angle += FULL_CIRC();
    }
    angle
}

/// hav(θ) = sin²(θ/2) = (1 - cos(θ)) / 2
pub(crate) fn haversine(theta: Angle) -> f64 {
    (1.0 - theta.get::<radian>().cos()) / 2.0
}

/// hav(θ) = sin²(θ/2) = (1 - cos(θ)) / 2
/// <=> cos(θ) = 1 - 2 * hav(θ)
/// => archav(x) = arccos(1 - 2x)
pub(crate) fn arcushaversine(h: f64) -> Angle {
    safe_acos(1.0 - 2.0 * h)
}

/// Saves acos from being called with an argument outside of its definition range due to numerical instability
pub(crate) fn safe_acos(cosine_argument: f64) -> Angle {
    if cosine_argument > 1. {
        return ANGLE_ZERO();
    } else if cosine_argument < -1. {
        return HALF_CIRC();
    }
    Angle::new::<radian>(cosine_argument.acos())
}

#[cfg(test)]
pub(crate) mod test {
    use super::*;

    pub(crate) fn THREE_QUARTER_CIRC() -> Angle {
        Angle::new::<radian>(2. * PI * 3. / 4.)
    }
    pub(crate) fn ONE_THIRD_CIRC() -> Angle {
        Angle::new::<radian>(2. * PI / 3.)
    }
    pub(crate) fn TWO_THIRDS_CIRC() -> Angle {
        Angle::new::<radian>(2. * PI * 2. / 3.)
    }
    pub(crate) fn ARCSEC() -> Angle {
        Angle::new::<radian>(2. * PI / (360. * 60. * 60.))
    }
    pub(crate) fn SECOND_ANGLE() -> Angle {
        Angle::new::<radian>(2. * PI / (24. * 60. * 60.))
    }

    pub(crate) fn angle_from_arcsecs(arcsec: f64) -> Angle {
        arcsec * ARCSEC()
    }

    pub(crate) fn angle_to_arcsecs(angle: &Angle) -> f64 {
        angle / &ARCSEC()
    }

    pub(crate) fn angle_from_second_angle(second_angle: f64) -> Angle {
        second_angle * SECOND_ANGLE()
    }
}
