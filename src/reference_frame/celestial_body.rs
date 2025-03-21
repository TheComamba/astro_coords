use std::fmt::Display;

use uom::{
    fmt::DisplayStyle,
    si::{
        angle::{degree, radian},
        angular_velocity::radian_per_second,
        f64::{Angle, AngularVelocity, Time},
        time::day,
    },
};

use crate::{
    angle_helper::{full_circ, normalized_angle},
    spherical::Spherical,
};

/// This enum contains the celestial bodies for which the equatorial reference frame can be defined.
///
/// The elements can be converted to the right ascension and declination of the north pole of the celestial body.
///
/// It is possible to define a custom celestial body by providing the right ascension and declination of its north pole, as well as its prime meridian offset at the epoch and its siderial angular velocity.
///
/// # Examples
/// ```
/// use uom::si::f64::Angle;
/// use uom::si::f64::Time;
/// use uom::si::angle::degree;
/// use uom::si::time::day;
/// use astro_coords::reference_frame::{CelestialBody, RotationalElements};
/// use astro_coords::spherical::Spherical;
///
/// let earth = CelestialBody::Earth;
/// let z = earth.z_axis();
/// let (ra, dec) = (z.longitude, z.latitude);
/// assert!(ra.get::<degree>() < 1e-5);
/// assert!(dec.get::<degree>() - 90.0 < 1e-5);
///
/// let custom_ra = Angle::new::<degree>(27.5);
/// let custom_dec = Angle::new::<degree>(119.24);
/// let custom_z = Spherical::new(custom_ra, custom_dec);
/// let custom_prime_meridian_offset = Angle::new::<degree>(23.0);
/// let custom_prime_meridian_rate = Angle::new::<degree>(40.0) / Time::new::<day>(1.0);
/// let custom_rotational_elements = RotationalElements {
///     z_axis: custom_z,
///     prime_meridian_offset_offset: custom_prime_meridian_offset,
///     prime_meridian_offset_rate: custom_prime_meridian_rate.into(),
/// };
/// let custom_body = CelestialBody::Custom(custom_rotational_elements);
/// let z = custom_body.z_axis();
/// let (ra, dec) = (z.longitude, z.latitude);
/// assert!((ra-custom_ra).get::<degree>() < 1e-5);
/// assert!((dec-custom_dec).get::<degree>() < 1e-5);
/// let time = Time::new::<day>(1.0);
/// let offset = custom_body.prime_meridian_offset(time);
/// let expected_offset = Angle::new::<degree>(23.0) + Angle::new::<degree>(40.0);
/// assert!((offset-expected_offset).get::<degree>().abs() < 1e-5);
/// ```
#[derive(Clone, Copy, Debug)]
pub enum CelestialBody {
    /// A celestial body with an arbitrary north-pole, provided as Right Ascension and Declination in Earth-Equatorial coordinates.
    Custom(RotationalElements),
    /// The central body in our solar system.
    Sun,
    /// The innermost and smallest planet in our solar system.
    Mercury,
    /// The second system in our solar system.
    Venus,
    /// The third planet in our solar system, and probably where you are right now.
    Earth,
    /// The fourth and outermost known rocky planet in our solar system.
    Mars,
    /// The fifth and largest planet in our solar system.
    Jupiter,
    /// The sixth planet in our solar system, known for its prominent rings.
    Saturn,
    /// The sevenths planet in our solar system, a gaseous ice-giant.
    Uranus,
    /// The eighths and outermost known planet in our solar system.
    Neptune,
}

#[derive(Clone, Copy, Debug)]
pub struct RotationalElements {
    /// The rotational axis, provided as Right Ascension and Declination in Earth-Equatorial coordinates.
    pub z_axis: Spherical,
    /// The prime meridian offset at the J2000 epoch.
    pub prime_meridian_offset_offset: Angle,
    /// The rate of change of the prime meridian offset.
    pub prime_meridian_offset_rate: AngularVelocity,
}

impl CelestialBody {
    /// Returns the right ascension and declination (in earth-equatorial coordinates) of the north pole of the celestial body.
    ///
    /// The data is taken from the [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf). Any time dependency is ignored.
    pub fn z_axis(&self) -> Spherical {
        match self {
            CelestialBody::Custom(rotational_elements) => rotational_elements.z_axis,
            CelestialBody::Sun => {
                Spherical::new(Angle::new::<degree>(286.13), Angle::new::<degree>(63.87))
            }
            CelestialBody::Mercury => Spherical::new(
                Angle::new::<degree>(281.0103),
                Angle::new::<degree>(61.4155),
            ),
            CelestialBody::Venus => {
                Spherical::new(Angle::new::<degree>(272.76), Angle::new::<degree>(67.16))
            }
            CelestialBody::Earth => {
                Spherical::new(Angle::new::<degree>(0.0), Angle::new::<degree>(90.0))
            }
            CelestialBody::Mars => Spherical::new(
                Angle::new::<degree>(317.269202),
                Angle::new::<degree>(54.432516),
            ),
            CelestialBody::Jupiter => Spherical::new(
                Angle::new::<degree>(268.056595),
                Angle::new::<degree>(64.495303),
            ),
            CelestialBody::Saturn => {
                Spherical::new(Angle::new::<degree>(40.589), Angle::new::<degree>(83.537))
            }
            CelestialBody::Uranus => {
                Spherical::new(Angle::new::<degree>(257.311), Angle::new::<degree>(-15.175))
            }
            CelestialBody::Neptune => {
                Spherical::new(Angle::new::<degree>(299.36), Angle::new::<degree>(43.46))
            }
        }
    }

    /// Returns the prime meridian offset of the reference frame at a given time since epoch.
    ///
    /// The prime meridian offset is the longitudal angle offset of the x-axis of the reference frame, relative to the point Q that the equatorial x-axis points at after rotating to the new z-axis. It is the current longitudial offset from the direction of earth's vernal equinox, or the x-direction in equatorial coordinates. It is denoted W in [Fig. 1 of the IAU Report](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf).
    /// The data here is taken from the [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf), except for that of earth, which is taken from Eq. 5.14 of [chapter 5 of technical note 36 of the international earth rotation and reference systems service](https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36_043.pdf?__blob=publicationFile&v=1).
    /// Only the constant and linear part of the time dependency is considered.
    /// The J2000 epoch for celestial bodies is defined as the 1st of January 2000, 12:00:00 TT.
    ///
    /// # Example
    /// ```
    /// use uom::si::f64::Time;
    /// use uom::si::angle::degree;
    /// use uom::si::time::day;
    /// use astro_coords::reference_frame::CelestialBody;
    ///
    /// // Vernal equinox was on the 20th of March 2000, which is 79 days after the J2000 epoch.
    /// let time_between_j2000_and_vernal_equinox = Time::new::<day>(79.);
    /// // At the vernal equinox, the prime meridian offset of the earth is roughly 0.
    /// let prime_meridian_offset = CelestialBody::Earth.prime_meridian_offset(time_between_j2000_and_vernal_equinox);
    /// // It is not exactly 0, because it is not guaranteed that the moment of vernal equinox is exactly at 12:00:00 TT.
    /// let acceptable_error_in_degrees = 2.;
    /// assert!(prime_meridian_offset.get::<degree>().abs() < acceptable_error_in_degrees);
    /// ```
    pub fn prime_meridian_offset(&self, time_since_epoch: Time) -> Angle {
        let offset = match self {
            CelestialBody::Custom(rotational_elements) => {
                rotational_elements.prime_meridian_offset_offset
                    + (rotational_elements.prime_meridian_offset_rate * time_since_epoch).into()
            }
            CelestialBody::Sun => {
                Angle::new::<degree>(84.176)
                    + Angle::new::<degree>(14.1844000) * time_since_epoch.get::<day>()
            }
            CelestialBody::Mercury => {
                Angle::new::<degree>(329.5988)
                    + Angle::new::<degree>(6.1385108) * time_since_epoch.get::<day>()
            }
            CelestialBody::Venus => {
                Angle::new::<degree>(160.20)
                    - Angle::new::<degree>(1.4813688) * time_since_epoch.get::<day>()
            }
            CelestialBody::Earth => {
                full_circ()
                    * (0.7790572732640 + 1.00273781191135448 * time_since_epoch.get::<day>())
            }
            CelestialBody::Mars => {
                Angle::new::<degree>(176.049863)
                    + Angle::new::<degree>(350.891982443297) * time_since_epoch.get::<day>()
            }
            CelestialBody::Jupiter => {
                Angle::new::<degree>(284.95)
                    + Angle::new::<degree>(870.5360000) * time_since_epoch.get::<day>()
            }
            CelestialBody::Saturn => {
                Angle::new::<degree>(38.90)
                    + Angle::new::<degree>(810.7939024) * time_since_epoch.get::<day>()
            }
            CelestialBody::Uranus => {
                Angle::new::<degree>(203.81)
                    - Angle::new::<degree>(501.1600928) * time_since_epoch.get::<day>()
            }
            CelestialBody::Neptune => {
                Angle::new::<degree>(249.978)
                    + Angle::new::<degree>(541.1397757) * time_since_epoch.get::<day>()
            }
        };
        normalized_angle(offset)
    }
}

impl Display for CelestialBody {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CelestialBody::Custom(rot) => {
                write!(f, "Custom body rotating around the {} axis, with a prime meridian offset of {}+{}*t",
                rot.z_axis,
                rot.prime_meridian_offset_offset.into_format_args(radian, DisplayStyle::Abbreviation),
                rot.prime_meridian_offset_rate.into_format_args(radian_per_second, DisplayStyle::Abbreviation))
            }
            CelestialBody::Sun => write!(f, "Sun"),
            CelestialBody::Mercury => write!(f, "Mercury"),
            CelestialBody::Venus => write!(f, "Venus"),
            CelestialBody::Earth => write!(f, "Earth"),
            CelestialBody::Mars => write!(f, "Mars"),
            CelestialBody::Jupiter => write!(f, "Jupiter"),
            CelestialBody::Saturn => write!(f, "Saturn"),
            CelestialBody::Uranus => write!(f, "Uranus"),
            CelestialBody::Neptune => write!(f, "Neptune"),
        }
    }
}
