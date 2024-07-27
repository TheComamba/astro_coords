use simple_si_units::geometry::Angle;
use std::f64::consts::PI;

pub(crate) const ANGLE_ZERO: Angle<f64> = Angle { rad: 0. };
pub(crate) const FULL_CIRC: Angle<f64> = Angle { rad: 2. * PI };
pub(crate) const QUARTER_CIRC: Angle<f64> = Angle { rad: 2. * PI / 4. };
pub(crate) const HALF_CIRC: Angle<f64> = Angle { rad: 2. * PI / 2. };

pub(crate) const DEGREE: Angle<f64> = Angle {
    rad: 2. * PI / 360.,
};
pub(crate) const EARTH_AXIS_TILT: Angle<f64> = Angle {
    rad: 23.439_281 * DEGREE.rad,
};

pub(crate) fn angle_eq_within(
    actual: Angle<f64>,
    expected: Angle<f64>,
    accuracy: Angle<f64>,
) -> bool {
    let diff = normalized_angle(actual - expected);
    diff.rad.abs() < accuracy.rad
}

/// Normalize the angle to a range of −π to +π radians, -180° to 180°.
pub fn normalized_angle(mut angle: Angle<f64>) -> Angle<f64> {
    angle.rad %= FULL_CIRC.rad;
    if angle > HALF_CIRC {
        angle -= FULL_CIRC;
    } else if angle < -HALF_CIRC {
        angle += FULL_CIRC;
    }
    angle
}

/// hav(θ) = sin²(θ/2) = (1 - cos(θ)) / 2
pub(crate) fn haversine(theta: Angle<f64>) -> f64 {
    (1.0 - theta.rad.cos()) / 2.0
}

/// hav(θ) = sin²(θ/2) = (1 - cos(θ)) / 2
/// <=> cos(θ) = 1 - 2 * hav(θ)
/// => archav(x) = arccos(1 - 2x)
pub(crate) fn arcushaversine(h: f64) -> Angle<f64> {
    safe_acos(1.0 - 2.0 * h)
}

/// Saves acos from being called with an argument outside of its definition range due to numerical instability
pub(crate) fn safe_acos(cosine_argument: f64) -> Angle<f64> {
    if cosine_argument > 1. {
        return ANGLE_ZERO;
    } else if cosine_argument < -1. {
        return HALF_CIRC;
    }
    Angle::from_radians(cosine_argument.acos())
}

#[cfg(test)]
pub(crate) mod test {
    use super::*;

    pub(crate) const THREE_QUARTER_CIRC: Angle<f64> = Angle {
        rad: 2. * PI * 3. / 4.,
    };
    pub(crate) const ONE_THIRD_CIRC: Angle<f64> = Angle { rad: 2. * PI / 3. };
    pub(crate) const TWO_THIRDS_CIRC: Angle<f64> = Angle {
        rad: 2. * PI * 2. / 3.,
    };
    pub(crate) const ARCSEC: Angle<f64> = Angle {
        rad: 2. * PI / (360. * 60. * 60.),
    };
    pub(crate) const SECOND_ANGLE: Angle<f64> = Angle {
        rad: 2. * PI / (24. * 60. * 60.),
    };

    pub(crate) fn angle_from_arcsecs(arcsec: f64) -> Angle<f64> {
        arcsec * ARCSEC
    }

    pub(crate) fn angle_to_arcsecs(angle: &Angle<f64>) -> f64 {
        angle / &ARCSEC
    }

    pub(crate) fn angle_from_second_angle(second_angle: f64) -> Angle<f64> {
        second_angle * SECOND_ANGLE
    }
}
