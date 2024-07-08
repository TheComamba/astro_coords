//! Right Ascension and Declination types and conversions.

use simple_si_units::geometry::Angle;
use std::fmt::Display;

pub struct RightAscension {
    pub hours: u8,
    pub minutes: u8,
    pub seconds: f64,
}

pub struct Declination {
    pub sign: Sgn,
    pub degrees: u8,
    pub minutes: u8,
    pub seconds: f64,
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
    pub const fn new(hours: u8, minutes: u8, seconds: f64) -> Self {
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

impl Display for RightAscension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:02}h{:02}m{:02}s",
            self.hours, self.minutes, self.seconds
        )
    }
}

impl Declination {
    pub const fn new(sign: Sgn, degrees: u8, minutes: u8, seconds: f64) -> Self {
        Self {
            sign,
            degrees,
            minutes,
            seconds,
        }
    }

    pub fn to_angle(&self) -> Angle<f64> {
        let sign = match self.sign {
            Sgn::Pos => 1.,
            Sgn::Neg => -1.,
        };
        let degrees = self.degrees as f64;
        let minutes = self.minutes as f64;
        let seconds = self.seconds as f64;

        sign * Angle::from_degrees(degrees + minutes / 60. + seconds / 3600.)
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
            "{}{:02}°{:02}'{:02}\"",
            sign, self.degrees, self.minutes, self.seconds
        )
    }
}

#[cfg(test)]
mod tests {
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
            Angle { rad: 1e-5 }
        ));
    }

    #[test]
    fn dec_one_second() {
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
