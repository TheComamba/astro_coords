use std::f64::consts::PI;

use uom::si::{
    angle::{degree, radian},
    f64::Angle,
};

#[inline]
pub(crate) fn angle_zero() -> Angle {
    Angle::new::<radian>(0.)
}
#[inline]
pub(crate) fn full_circ() -> Angle {
    Angle::new::<radian>(2. * PI)
}
#[inline]
pub(crate) fn quarter_circ() -> Angle {
    Angle::new::<radian>(2. * PI / 4.)
}
#[inline]
pub(crate) fn half_circ() -> Angle {
    Angle::new::<radian>(2. * PI / 2.)
}
#[inline]
pub(crate) fn earth_axis_tilt() -> Angle {
    Angle::new::<degree>(23.439_281)
}
#[inline]
pub(crate) fn galactic_longitude_of_north_celestial_pole() -> Angle {
    Angle::new::<radian>(122.93314)
}
#[inline]
pub(crate) fn ra_of_galactic_north() -> Angle {
    Angle::new::<radian>(3.3658674624710652)
}
#[inline]
pub(crate) fn dec_of_galactic_north() -> Angle {
    Angle::new::<radian>(27.13)
}

pub(crate) fn angle_eq_within(actual: Angle, expected: Angle, accuracy: Angle) -> bool {
    let diff = normalized_angle(actual - expected);
    diff.get::<radian>().abs() < accuracy.get::<radian>()
}

/// Normalize the angle to a range of −π to +π radians, -180° to 180°.
pub fn normalized_angle(mut angle: Angle) -> Angle {
    angle %= full_circ();
    if angle > half_circ() {
        angle -= full_circ();
    } else if angle < -half_circ() {
        angle += full_circ();
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
        return angle_zero();
    } else if cosine_argument < -1. {
        return half_circ();
    }
    Angle::new::<radian>(cosine_argument.acos())
}

#[cfg(test)]
pub(crate) mod test {
    use super::*;

    #[inline]
    pub(crate) fn three_quarter_circ() -> Angle {
        Angle::new::<radian>(2. * PI * 3. / 4.)
    }
    #[inline]
    pub(crate) fn one_third_circ() -> Angle {
        Angle::new::<radian>(2. * PI / 3.)
    }
    #[inline]
    pub(crate) fn two_thirds_circ() -> Angle {
        Angle::new::<radian>(2. * PI * 2. / 3.)
    }
    #[inline]
    pub(crate) fn arcsec() -> Angle {
        Angle::new::<radian>(2. * PI / (360. * 60. * 60.))
    }
    #[inline]
    pub(crate) fn second_angle() -> Angle {
        Angle::new::<radian>(2. * PI / (24. * 60. * 60.))
    }
    #[inline]
    pub(crate) fn angle_from_arcsecs(arcsecs: f64) -> Angle {
        arcsecs * arcsec()
    }
    #[inline]
    pub(crate) fn angle_to_arcsecs(angle: &Angle) -> f64 {
        angle.get::<radian>() / arcsec().get::<radian>()
    }
    #[inline]
    pub(crate) fn angle_from_second_angle(second_angles: f64) -> Angle {
        second_angles * second_angle()
    }
}
