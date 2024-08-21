//! This module contains the `ReferenceFrame` enum, which connects mathematical coordinates to physical positions and directions.

use std::fmt::Display;

use simple_si_units::geometry::Angle;

use crate::{
    angle_helper::{EARTH_AXIS_TILT, QUARTER_CIRC},
    ra_and_dec::RightAscension,
    spherical::Spherical,
};

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

impl ReferenceFrame {
    /// Returns the z-axis, or the direction orthogonal to the fundamental plane, of the reference frame.
    ///
    /// The z-axis is expessed as a spherical direction in the equatorial frame of reference.
    ///
    /// # Example
    /// ```
    /// use astro_coords::reference_frame::ReferenceFrame;
    /// use astro_coords::spherical::Spherical;
    /// use astro_coords::ra_and_dec::RightAscension;
    /// use astro_coords::traits::*;
    /// use simple_si_units::geometry::Angle;
    ///
    /// // The z-axis of the equatorial frame expressed in equatorial coordinates is, well, the z-direction.
    /// let equatorial_north = Spherical::Z_DIRECTION;
    ///
    /// // The x-axis in both equatorial and ecliptic coordinates is the direction of the vernal equinox.
    /// // The y-axis on the other hand it rotated around the x-axis by the axial tilt of the Earth.
    /// let earth_axial_tilt = Angle::from_degrees(23.439);
    /// let ecliptic_north = Spherical::Z_DIRECTION.rotated_x(earth_axial_tilt);
    /// // After a quarter-year, the sun is high in the sky on the northern hemisphere.
    /// // Because it encircles the earth in mathematically positive direction, this means that the y-component of ecliptic_north is negative.
    /// assert!(ecliptic_north.to_direction().y() < 0., "{}", ecliptic_north);
    ///
    /// // RA and Dec of the north pole of the Milky Way galaxy can be found online.
    /// let galactic_north_right_ascension = RightAscension::new(12, 49, 0.).to_angle();
    /// let galactic_north_declination = Angle::from_degrees(27.4);
    /// let galactic_north = Spherical::new(galactic_north_right_ascension, galactic_north_declination);
    ///
    /// let acc = Angle::from_degrees(1e-5);
    /// assert!(ReferenceFrame::Equatorial.z_axis().eq_within(&equatorial_north, acc), "{}", ReferenceFrame::Equatorial.z_axis());
    /// assert!(ReferenceFrame::Ecliptic.z_axis().eq_within(&ecliptic_north, acc), "{}", ReferenceFrame::Ecliptic.z_axis());
    /// assert!(ReferenceFrame::Galactic.z_axis().eq_within(&galactic_north, acc), "{}", ReferenceFrame::Galactic.z_axis());
    /// ```
    pub fn z_axis(&self) -> Spherical {
        match self {
            ReferenceFrame::Equatorial => Spherical::Z_DIRECTION,
            ReferenceFrame::Ecliptic => {
                Spherical::new(-QUARTER_CIRC, QUARTER_CIRC - EARTH_AXIS_TILT)
            }
            ReferenceFrame::Galactic => Spherical::new(
                RightAscension::new(12, 49, 0.).to_angle(),
                Angle::from_degrees(27.4),
            ),
        }
    }

    /// Returns the prime meridian offset of the reference frame.
    ///
    /// The prime meridian offset is the longitudal angle offset of the x-axis of the reference frame, relative to the point Q that the equatorial x-axis points at after rotating to the new z-axis. It is denoted W in [Fig. 1 of the IAU Report](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf).
    ///
    /// # Example
    /// ```
    /// use astro_coords::reference_frame::ReferenceFrame;
    /// use astro_coords::direction::Direction;
    /// use astro_coords::spherical::Spherical;
    /// use astro_coords::ra_and_dec::RightAscension;
    /// use astro_coords::traits::*;
    /// use simple_si_units::geometry::Angle;
    ///
    /// // The equatorial reference frame has no offset to itself.
    /// assert!((ReferenceFrame::Equatorial.prime_meridian_offset() - Angle::from_deg(0.)).rad.abs() < 1e-5);
    ///
    /// // The ecliptic x-axis points in the same direction (the vernal equinox) as the equatorial x-axis.
    /// assert!((ReferenceFrame::Ecliptic.prime_meridian_offset() - Angle::from_deg(0.)).rad.abs() < 1e-5);
    ///
    /// // Galactic north in equatorial coordinates is defined as a=12h49m, d=27.4°.
    /// let galactic_north_ra = RightAscension::new(12, 49, 0.).to_angle();
    /// // The x-axis of the galactic reference frame points towards the galactic center, which for the purpose of the reference frame is fixed at a=17h42m24s, d=-28.92°.
    /// let galactic_center_ra = RightAscension::new(17, 42, 24.).to_angle();
    /// let galactic_center_dec = Angle::from_deg(-28.92);
    /// let galactic_center = Spherical::new(galactic_center_ra, galactic_center_dec).to_direction();
    /// // The point Q is found by rotating the x-axis by 90°-12h49m (compare [Fig. 1 of the IAU Report](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf)).
    /// let Q = Direction::X.rotated_z(Angle::from_deg(90.) - galactic_north_ra);
    /// // The prime meridian offset is the angle between the galactic center and the point Q.
    /// let expected_offset = galactic_center.angle_to(&Q);
    /// assert!((ReferenceFrame::Galactic.prime_meridian_offset() - expected_offset).rad.abs() < 1e-5);
    /// ```
    pub fn prime_meridian_offset(&self) -> Angle<f64> {
        match self {
            ReferenceFrame::Equatorial => Angle::from_deg(0.),
            ReferenceFrame::Ecliptic => Angle::from_deg(0.),
            ReferenceFrame::Galactic => Angle::from_deg(0.),
        }
    }
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
