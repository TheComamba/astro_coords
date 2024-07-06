use simple_si_units::geometry::Angle;
use std::f64::consts::PI;

pub(crate) const ANGLE_ZERO: Angle<f64> = Angle { rad: 0. };
pub(crate) const FULL_CIRC: Angle<f64> = Angle { rad: 2. * PI };
pub(crate) const QUARTER_CIRC: Angle<f64> = Angle { rad: 2. * PI / 4. };
pub(crate) const HALF_CIRC: Angle<f64> = Angle { rad: 2. * PI / 2. };
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
}

pub(crate) const DEGREE: Angle<f64> = Angle {
    rad: 2. * PI / 360.,
};
pub(crate) const ARCSEC: Angle<f64> = Angle {
    rad: 2. * PI / (360. * 60. * 60.),
};
pub(crate) const SECOND_ANGLE: Angle<f64> = Angle {
    rad: 2. * PI / (24. * 60. * 60.),
};
pub(crate) const EARTH_AXIS_TILT: Angle<f64> = Angle {
    rad: 23.439_281 * DEGREE.rad,
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

pub(crate) fn angle_to_second_angle(angle: &Angle<f64>) -> f64 {
    angle / &SECOND_ANGLE
}

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
