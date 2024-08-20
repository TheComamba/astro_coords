//! This module contains the `ReferenceFrame` enum, which connects mathematical coordinates to physical positions and directions.

use std::fmt::Display;

use simple_si_units::geometry::Angle;

use crate::{
    angle_helper::{EARTH_AXIS_TILT, QUARTER_CIRC},
    direction::Direction,
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
    /// The z-axis is expessed as a direction in the equatorial frame of reference.
    ///
    /// # Example
    /// ```
    /// use astro_coords::reference_frame::ReferenceFrame;
    /// use astro_coords::direction::Direction;
    /// use astro_coords::spherical::Spherical;
    /// use astro_coords::ra_and_dec::RightAscension;
    /// use simple_si_units::angle::Angle;
    ///
    /// // The z-axis of the equatorial frame expressed in equatorial coordinates is, well, the z-direction.
    /// let equatorial_north = Direction::Z;
    ///
    /// // The x-axis in both equatorial and ecliptic coordinates is the direction of the vernal equinox.
    /// // The y-axis on the other hand it rotated around the x-axis by the axial tilt of the Earth.
    /// let earth_axial_tilt = Angle::from_degrees(23.439);
    /// let ecliptic_north = Direction::Z.rotated_x(earth_axial_tilt);
    /// // After a quarter-year, the sun is high in the sky on the northern hemisphere.
    /// // Because it encircles the earth in mathematically positive direction, this means that the y-component of ecliptic_north is negative.
    /// assert!(ecliptic_north.y() < 0.);
    ///
    /// // RA and Dec of the north pole of the Milky Way galaxy can be found online.
    /// let galactic_north_right_ascension = RightAscension::new(12, 51, 24.).to_angle();
    /// let galactic_north_declination = Angle::from_degrees(27.13);
    /// let galactic_north = Spherical::new(galactic_north_right_ascension, galactic_north_declination);
    /// let galactic_north = galactic_north.to_direction();
    ///
    /// assert_eq!(ReferenceFrame::Equatorial.z_axis().eq_within(equatorial_north, 1e-5));
    /// assert_eq!(ReferenceFrame::Ecliptic.z_axis().eq_within(ecliptic_north, 1e-5));
    /// assert_eq!(ReferenceFrame::Galactic.z_axis().eq_within(galactic_north, 1e-5));
    /// ```
    pub fn z_axis(&self) -> Direction {
        match self {
            ReferenceFrame::Equatorial => Direction::Z,
            ReferenceFrame::Ecliptic => {
                Spherical::new(QUARTER_CIRC, QUARTER_CIRC - EARTH_AXIS_TILT).to_direction()
            }
            ReferenceFrame::Galactic => Spherical::new(
                RightAscension::new(12, 51, 24.).to_angle(),
                Angle::from_degrees(27.13),
            )
            .to_direction(),
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
