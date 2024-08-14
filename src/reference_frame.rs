//! This module contains the `ReferenceFrame` enum, which connects mathematical coordinates to physical positions and directions.

use simple_si_units::geometry::Angle;

/// A reference frame provides the connection between mathematical coordinates (say, (0,0,1) or the z-direction) and physical positions or directions (say, the direction of the North Pole of earth).
///
/// This enum needs to be provided  to eliminate ambiguity when converting between different coordinate systems.
pub enum ReferenceFrame {
    Equatorial(CelestialBody),
    Ecliptic,
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
pub enum CelestialBody {
    Custom(Angle<f64>, Angle<f64>),
    Earth,
}

impl CelestialBody {
    /// Returns the right ascension and declination (in earth-equatorial coordinates) of the north pole of the celestial body.
    ///
    /// The data is taken from the [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf).
    pub fn get_ra_and_dec(&self) -> (Angle<f64>, Angle<f64>) {
        match self {
            CelestialBody::Custom(ra, dec) => (*ra, *dec),
            CelestialBody::Earth => (Angle::from_degrees(0.0), Angle::from_degrees(90.0)),
        }
    }
}
