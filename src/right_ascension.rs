use std::fmt::Display;

use simple_si_units::geometry::Angle;

use crate::astro_display::AstroDisplay;

pub struct RightAscension {
    pub(super) hours: i8,
    pub(super) minutes: i8,
    pub(super) seconds: i8,
}

impl RightAscension {
    pub const fn new(hours: i8, minutes: i8, seconds: i8) -> Self {
        Self {
            hours,
            minutes,
            seconds,
        }
    }

    pub fn to_angle(&self) -> Angle<f64> {
        let hours = self.hours as f64;
        let minutes = self.minutes as f64;
        let seconds = self.seconds as f64;

        Angle::from_degrees((hours + minutes / 60. + seconds / 3600.) * 15.)
    }
}

impl AstroDisplay for RightAscension {
    fn astro_display(&self) -> String {
        format!("{:02}h{:02}m{:02}s", self.hours, self.minutes, self.seconds)
    }
}

impl Display for RightAscension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.astro_display())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::angle::{angle_eq, angle_from_second_angle, angle_to_arcsecs};

    #[test]
    fn one_second() {
        let dec = RightAscension::new(0, 0, 1);
        let expected = angle_from_second_angle(1.);
        println!("{}", angle_to_arcsecs(&dec.to_angle()));
        println!("{}", angle_to_arcsecs(&expected));
        assert!(angle_eq(dec.to_angle(), expected));
    }
}
