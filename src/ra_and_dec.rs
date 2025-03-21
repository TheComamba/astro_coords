//! Right Ascension and Declination types and conversions.

use std::fmt::Display;
use uom::si::{angle::degree, f64::Angle};

/// Right ascension is the angular distance of a point eastward along the celestial equator from the vernal equinox to the point in question.
///
/// It is measured in hours, minutes, and seconds, with 24 hours corresponding to a full circle.
///
/// The main functionality of this struct is to convert from hours, minutes, and seconds to an angle.
///
/// # Examples
/// ```
/// use uom::si::angle::degree;
/// use astro_coords::ra_and_dec::RightAscension;
///
/// let ra = RightAscension::new(1, 2, 3.456);
/// assert_eq!(format!("{}", ra), "01h02m03.456s");
///
/// let ra = RightAscension::new(6, 0, 0.);
/// assert!((ra.to_angle().get::<degree>() - 90.).abs() < 1e-5);
///
/// let one_hour = RightAscension::new(1, 0, 0.);
/// let sixty_minutes = RightAscension::new(0, 60, 0.);
/// assert!((one_hour.to_angle() - sixty_minutes.to_angle()).get::<degree>().abs() < 1e-5);
///
/// let one_minute = RightAscension::new(0, 1, 0.);
/// let sixty_seconds = RightAscension::new(0, 0, 60.);
/// assert!((one_minute.to_angle() - sixty_seconds.to_angle()).get::<degree>().abs() < 1e-5);
/// ```
pub struct RightAscension {
    /// Subdivision of the equatorial plane into 24 hours.
    pub hours: u8,
    /// Subdivision of an hour into 60 minutes.
    pub minutes: u8,
    /// Subdivision of a minute into 60 seconds.
    pub seconds: f64,
}

/// Declination is the angular distance of a point north or south of the celestial equator.
///
/// It is measured in degrees, arcminutes, and arcseconds.
/// 90 degrees correspond to the north celestial pole and -90 degrees to the south celestial pole.
/// One arcminute is 1/60 of a degree, and one arcsecond is 1/60 of an arcminute.
///
/// The main functionality of this struct is to convert from degrees, arcminutes, and arcseconds to an angle.
///
/// # Examples
/// ```
/// use uom::si::angle::degree;
/// use astro_coords::ra_and_dec::{Declination, Sgn};
///
/// let dec = Declination::new(Sgn::Pos, 1, 2, 3.456);
/// assert_eq!(format!("{}", dec), "+01°02'03.456\"");
///
/// let dec = Declination::new(Sgn::Neg, 45, 0, 0.);
/// assert!((dec.to_angle().get::<degree>() + 45.).abs() < 1e-5);
///
/// let one_degree = Declination::new(Sgn::Pos, 1, 0, 0.);
/// let sixty_arcminutes = Declination::new(Sgn::Pos, 0, 60, 0.);
/// assert!((one_degree.to_angle() - sixty_arcminutes.to_angle()).get::<degree>().abs() < 1e-5);
///
/// let one_arcminute = Declination::new(Sgn::Pos, 0, 1, 0.);
/// let sixty_arcseconds = Declination::new(Sgn::Pos, 0, 0, 60.);
/// assert!((one_arcminute.to_angle() - sixty_arcseconds.to_angle()).get::<degree>().abs() < 1e-5);
/// ```
pub struct Declination {
    /// Sign of the declination, denoting northern or southern hemisphere.
    pub sign: Sgn,
    /// The latitude of the point in degrees.
    pub degrees: u8,
    /// Subdivision of a degree into 60 arcminutes.
    pub arcminutes: u8,
    /// Subdivision of an arcminute into 60 arcseconds.
    pub arcseconds: f64,
}

/// Sign of a declination
///
/// This simple enum resolves the ambiguity that would arise from 0 == -0 in case degrees were stored in a signed integer.
#[derive(Copy, Clone)]
pub enum Sgn {
    /// Positive sign, corresponding to the northern hemisphere.
    Pos,
    /// Negative sign, corresponding to the southern hemisphere.
    Neg,
}

impl RightAscension {
    /// Create a new RightAscension instance.
    pub const fn new(hours: u8, minutes: u8, seconds: f64) -> Self {
        Self {
            hours,
            minutes,
            seconds,
        }
    }

    /// Convert the right ascension to an angle.
    pub fn to_angle(&self) -> Angle {
        let hours = self.hours as f64;
        let minutes = self.minutes as f64;
        let seconds = self.seconds as f64;

        Angle::new::<degree>((hours + minutes / 60. + seconds / 3600.) * 15.)
    }
}

impl Display for RightAscension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:02}h{:02}m{:06.3}s",
            self.hours, self.minutes, self.seconds
        )
    }
}

impl Declination {
    /// Create a new Declination instance.
    pub const fn new(sign: Sgn, degrees: u8, arcminutes: u8, arcseconds: f64) -> Self {
        Self {
            sign,
            degrees,
            arcminutes,
            arcseconds,
        }
    }

    /// Convert the declination to an angle.
    pub fn to_angle(&self) -> Angle {
        let sign = match self.sign {
            Sgn::Pos => 1.,
            Sgn::Neg => -1.,
        };
        let degrees = self.degrees as f64;
        let minutes = self.arcminutes as f64;
        let seconds = self.arcseconds as f64;

        sign * Angle::new::<degree>(degrees + minutes / 60. + seconds / 3600.)
    }
}

impl Display for Declination {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let sign = match self.sign {
            Sgn::Pos => "+",
            Sgn::Neg => "-",
        };
        write!(
            f,
            "{}{:02}°{:02}'{:06.3}\"",
            sign, self.degrees, self.arcminutes, self.arcseconds
        )
    }
}

#[cfg(test)]
mod tests {
    use uom::si::angle::radian;

    use crate::angle_helper::{angle_eq_within, test::*};

    use super::*;

    #[test]
    fn ra_one_second() {
        let dec = RightAscension::new(0, 0, 1.);
        let expected = angle_from_second_angle(1.);
        println!("{}", angle_to_arcsecs(&dec.to_angle()));
        println!("{}", angle_to_arcsecs(&expected));
        assert!(angle_eq_within(
            dec.to_angle(),
            expected,
            Angle::new::<radian>(1e-5)
        ));
    }

    #[test]
    fn dec_one_arcsecond() {
        let dec = Declination::new(Sgn::Pos, 0, 0, 1.);
        assert!((angle_to_arcsecs(&dec.to_angle()) - 1.) < 1e-5);
    }

    #[test]
    fn dec_small_angle_roundtrips() {
        const STEPS: u8 = 10;
        for sec in 0..STEPS {
            for min in 0..STEPS {
                for sign in [Sgn::Pos, Sgn::Neg] {
                    let dec = Declination::new(sign, 0, min, sec as f64);
                    let angle_abs = ((min as u32) * 60 + sec as u32) as f64;
                    let sign = match sign {
                        Sgn::Pos => 1.,
                        Sgn::Neg => -1.,
                    };
                    let angle = angle_from_arcsecs(sign * angle_abs);
                    println!("{} sign {} min {} sec", sign, min, sec);
                    println!("expected: \n{} sec", angle_to_arcsecs(&angle));
                    println!("actual: \n{} sec", angle_to_arcsecs(&dec.to_angle()));
                    let diff = (angle_to_arcsecs(&dec.to_angle()) - angle_to_arcsecs(&angle)).abs();
                    println!("diff: \n{} sec", diff);
                    let mut accuracy = 1e-5 * angle_abs;
                    if accuracy < 1e-5 {
                        accuracy = 1e-5;
                    }
                    println!("accuracy: \n{} sec", accuracy);
                    assert!(diff < accuracy);
                }
            }
        }
    }
}
