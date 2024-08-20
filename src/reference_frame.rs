//! This module contains the `ReferenceFrame` enum, which connects mathematical coordinates to physical positions and directions.

use std::fmt::Display;

use simple_si_units::geometry::Angle;

/// A reference frame provides the connection between mathematical coordinates (say, (0,0,1) or the z-direction) and physical positions or directions (say, the direction of the North Pole of earth).
///
/// This enum needs to be provided  to eliminate ambiguity when converting between different coordinate systems.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReferenceFrame {
    /// The Earth-Equatorial reference frame.
    ///
    /// Equatorial coordinates are oriented such that the z-axis points along the Earth's axis of rotation or north pole. The remaining orientation ambiguity is resolved by any of the four following facts:
    /// - At the vernal equinox, the x-axis points directly towards the sun.
    /// - At the summer solstice, when the sun appears high in the sky on the northern hemisphere, its y-projection and z-projection are both positive.
    /// - At the autumn equinox, the x-axis points directly away from the sun.
    /// - At the winter solstice, when the sun appears high in the sky on the southern hemisphere, its y-projection and z-projection are both negative.
    ///
    /// Note that the orientation of the coordinate system does not change with time, even though the Earth is rotating.
    Equatorial,
    /// The Earth-Ecliptic reference frame.
    ///
    /// Ecliptic coordinates are oriented such that the z-axis points along the normal of the plane of the Earth's orbit around the sun. The x-axis points towards the vernal equinox (the same as for Equatorial coordinates), and the y-axis completes the right-handed coordinate system.
    Ecliptic,
    /// The Galactic reference frame.
    ///
    /// Galactic coordinates are oriented such that the z-axis points along the normal of the plane of the Milky Way galaxy. The x-axis points towards the galactic center, and the y-axis completes the right-handed coordinate system.
    Galactic,
}

impl Display for ReferenceFrame {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReferenceFrame::Equatorial => write!(f, "Equatorial"),
            ReferenceFrame::Ecliptic => write!(f, "Ecliptic"),
            ReferenceFrame::Galactic => write!(f, "Galactic"),
        }
    }
}

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
/// let (ra, dec) = earth.get_ra_and_dec();
/// assert!(ra.to_degrees() < 1e-5);
/// assert!(dec.to_degrees() - 90.0 < 1e-5);
///
/// let custom_ra = Angle::from_degrees(27.5);
/// let custom_dec = Angle::from_degrees(119.24);
/// let custom_body = CelestialBody::Custom(custom_ra, custom_dec);
/// let (ra, dec) = custom_body.get_ra_and_dec();
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
    pub fn get_ra_and_dec(&self) -> (Angle<f64>, Angle<f64>) {
        match self {
            CelestialBody::Custom(ra, dec) => (*ra, *dec),
            CelestialBody::Sun => (Angle::from_deg(286.13), Angle::from_deg(63.87)),
            CelestialBody::Mercury => (Angle::from_deg(281.0103), Angle::from_deg(61.4155)),
            CelestialBody::Venus => (Angle::from_deg(272.76), Angle::from_deg(67.16)),
            CelestialBody::Earth => (Angle::from_deg(0.0), Angle::from_deg(90.0)),
            CelestialBody::Mars => (Angle::from_deg(317.269202), Angle::from_deg(54.432516)),
            CelestialBody::Jupiter => (Angle::from_deg(268.056595), Angle::from_deg(64.495303)),
            CelestialBody::Saturn => (Angle::from_deg(40.589), Angle::from_deg(83.537)),
            CelestialBody::Uranus => (Angle::from_deg(257.311), Angle::from_deg(-15.175)),
            CelestialBody::Neptune => (Angle::from_deg(299.36), Angle::from_deg(43.46)),
        }
    }
}

impl Display for CelestialBody {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CelestialBody::Custom(ra, dec) => write!(f, "Custom body with RA={} and Dec={}", ra, dec),
            CelestialBody::Sun => write!(f, "the Sun"),
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

impl PartialEq for CelestialBody {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (CelestialBody::Custom(_, _), _) | (_, CelestialBody::Custom(_, _)) => false,
            _ => std::mem::discriminant(self) == std::mem::discriminant(other),
        }
    }
}

impl Eq for CelestialBody {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn custom_celestial_bodies_are_never_equal() {
        let custom_ra = Angle::from_degrees(27.5);
        let custom_dec = Angle::from_degrees(119.24);
        let custom_body = CelestialBody::Custom(custom_ra, custom_dec);

        assert_ne!(custom_body, CelestialBody::Sun);
        assert_ne!(custom_body, CelestialBody::Custom(custom_ra, custom_dec));
    }

    #[test]
    fn non_custom_celestial_bodies_are_comparable() {
        assert_eq!(CelestialBody::Sun, CelestialBody::Sun);
        assert_eq!(CelestialBody::Earth, CelestialBody::Earth);
        assert_ne!(CelestialBody::Sun, CelestialBody::Earth);
    }
}
