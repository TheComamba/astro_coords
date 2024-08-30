use simple_si_units::{base::Time, geometry::Angle};
use std::fmt::Display;

use crate::{angle_helper::FULL_CIRC, spherical::Spherical};

/// This enum contains the celestial bodies for which the equatorial reference frame can be defined.
///
/// The elements can be converted to the right ascension and declination of the north pole of the celestial body.
///
/// It is possible to define a custom celestial body by providing the right ascension and declination of its north pole.
///
/// # Examples
/// ```
/// use simple_si_units::geometry::Angle;
/// use astro_coords::reference_frame::CelestialBody;
///
/// let earth = CelestialBody::Earth;
/// let z = earth.z_axis();
/// let (ra, dec) = (z.longitude, z.latitude);
/// assert!(ra.to_degrees() < 1e-5);
/// assert!(dec.to_degrees() - 90.0 < 1e-5);
///
/// let custom_ra = Angle::from_degrees(27.5);
/// let custom_dec = Angle::from_degrees(119.24);
/// let custom_body = CelestialBody::Custom(custom_ra, custom_dec);
/// let z = custom_body.z_axis();
/// let (ra, dec) = (z.longitude, z.latitude);
/// assert!((ra-custom_ra).to_degrees() < 1e-5);
/// assert!((dec-custom_dec).to_degrees() < 1e-5);
/// ```
#[derive(Clone, Copy, Debug)]
pub enum CelestialBody {
    /// A celestial body with an arbitrary north-pole, provided as Right Ascension and Declination in Earth-Equatorial coordinates.
    Custom(Angle<f64>, Angle<f64>),
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
    /// The eightths and outermost known planet in our solar system.
    Neptune,
}

impl CelestialBody {
    /// Returns the right ascension and declination (in earth-equatorial coordinates) of the north pole of the celestial body.
    ///
    /// The data is taken from the [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf). Any time dependency is ignored.
    pub fn z_axis(&self) -> Spherical {
        match self {
            CelestialBody::Custom(ra, dec) => Spherical::new(*ra, *dec),
            CelestialBody::Sun => Spherical::new(Angle::from_deg(286.13), Angle::from_deg(63.87)),
            CelestialBody::Mercury => {
                Spherical::new(Angle::from_deg(281.0103), Angle::from_deg(61.4155))
            }
            CelestialBody::Venus => Spherical::new(Angle::from_deg(272.76), Angle::from_deg(67.16)),
            CelestialBody::Earth => Spherical::new(Angle::from_deg(0.0), Angle::from_deg(90.0)),
            CelestialBody::Mars => {
                Spherical::new(Angle::from_deg(317.269202), Angle::from_deg(54.432516))
            }
            CelestialBody::Jupiter => {
                Spherical::new(Angle::from_deg(268.056595), Angle::from_deg(64.495303))
            }
            CelestialBody::Saturn => {
                Spherical::new(Angle::from_deg(40.589), Angle::from_deg(83.537))
            }
            CelestialBody::Uranus => {
                Spherical::new(Angle::from_deg(257.311), Angle::from_deg(-15.175))
            }
            CelestialBody::Neptune => {
                Spherical::new(Angle::from_deg(299.36), Angle::from_deg(43.46))
            }
        }
    }

    /// Returns the prime meridian offset of the reference frame.
    ///
    /// The prime meridian offset is the longitudal angle offset of the x-axis of the reference frame, relative to the point Q that the equatorial x-axis points at after rotating to the new z-axis. It is the current longitudial offset from the direction of earth's vernal equinox, or the x-direction in equatorial coordinates. It is denoted W in [Fig. 1 of the IAU Report](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf).
    /// The data here is taken from the [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf), except for that of earth, which is taken from Eq. 5.14 of [chapter 5 of technical note 36 of the international earth rotation and reference systems service](https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36_043.pdf?__blob=publicationFile&v=1).
    pub fn prime_meridian_offset(&self, time_since_epoch: Time<f64>) -> Angle<f64> {
        match self {
            CelestialBody::Sun => {
                Angle::from_deg(84.176) + Angle::from_deg(14.1844000) * time_since_epoch.to_days()
            }
            CelestialBody::Mercury => {
                Angle::from_deg(329.5988) + Angle::from_deg(6.1385108) * time_since_epoch.to_days()
            }
            CelestialBody::Venus => {
                Angle::from_deg(160.20) - Angle::from_deg(1.4813688) * time_since_epoch.to_days()
            }
            CelestialBody::Earth => {
                FULL_CIRC * (0.7790572732640 + 1.00273781191135448 * time_since_epoch.to_days())
            }
            CelestialBody::Mars => {
                Angle::from_deg(176.049863)
                    + Angle::from_deg(350.891982443297) * time_since_epoch.to_days()
            }
            CelestialBody::Jupiter => {
                Angle::from_deg(284.95) + Angle::from_deg(870.5360000) * time_since_epoch.to_days()
            }
            CelestialBody::Saturn => {
                Angle::from_deg(38.90) + Angle::from_deg(810.7939024) * time_since_epoch.to_days()
            }
            CelestialBody::Uranus => {
                Angle::from_deg(203.81) - Angle::from_deg(501.1600928) * time_since_epoch.to_days()
            }
            CelestialBody::Neptune => {
                Angle::from_deg(249.978) + Angle::from_deg(541.1397757) * time_since_epoch.to_days()
            }
            _ => todo!(),
        }
    }
}

impl Display for CelestialBody {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CelestialBody::Custom(ra, dec) => {
                write!(f, "Custom body with RA={} and Dec={}", ra, dec)
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
