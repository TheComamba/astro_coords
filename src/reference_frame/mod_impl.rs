//! This module contains the `ReferenceFrame` enum, which connects mathematical coordinates to physical positions and directions.

use std::fmt::Display;

use simple_si_units::geometry::Angle;

use crate::{
    angle_helper::{EARTH_AXIS_TILT, HALF_CIRC, QUARTER_CIRC},
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

impl Display for ReferenceFrame {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReferenceFrame::Equatorial => write!(f, "Equatorial"),
            ReferenceFrame::Ecliptic => write!(f, "Ecliptic"),
            ReferenceFrame::Galactic => write!(f, "Galactic"),
        }
    }
}
