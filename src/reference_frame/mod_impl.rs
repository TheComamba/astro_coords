//! This module contains the `ReferenceFrame` enum, which connects mathematical coordinates to physical positions and directions.

use std::fmt::Display;

use uom::si::f64::Time;

use super::CelestialBody;

/// A reference frame provides the connection between mathematical coordinates (say, (0,0,1) or the z-direction) and physical positions or directions (say, the direction of the North Pole of earth).
///
/// This enum needs to be provided  to eliminate ambiguity when converting between different coordinate systems.
#[derive(Clone, Copy, Debug)]
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
    /// The co-rotating reference frame for a specific celestial body at a specific time since the epoch.
    ///
    /// The z-axis points along the north pole of the celestial body, which is not necessarily to axis of positive rotation, but rather the axis of rotation that points towards the positive half-dome of the ecliptic.
    /// The x-axis points along a direction usually specified by tome distinguishable feature of the body, which is recommended by the IAU working group.
    Cartographic(CelestialBody, Time),
}

impl Display for ReferenceFrame {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReferenceFrame::Equatorial => write!(f, "Equatorial"),
            ReferenceFrame::Ecliptic => write!(f, "Ecliptic"),
            ReferenceFrame::Galactic => write!(f, "Galactic"),
            ReferenceFrame::Cartographic(body, _) => write!(f, "{} Cartographic", body),
        }
    }
}

impl PartialEq for ReferenceFrame {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (ReferenceFrame::Cartographic(_, _), _) | (_, ReferenceFrame::Cartographic(_, _)) => {
                false
            }
            _ => std::mem::discriminant(self) == std::mem::discriminant(other),
        }
    }
}

impl Eq for ReferenceFrame {}

#[cfg(test)]
mod tests {
    use uom::si::time::second;

    use super::*;

    #[test]
    fn cartographic_reference_frames_are_never_equal() {
        let body = CelestialBody::Sun;
        let time = Time::new::<second>(0.0);
        let cartographic_ref = ReferenceFrame::Cartographic(body, time);

        assert_ne!(cartographic_ref, ReferenceFrame::Equatorial);
        assert_ne!(cartographic_ref, cartographic_ref);
    }

    #[test]
    fn non_cartographic_reference_frames_are_comparable() {
        assert_eq!(ReferenceFrame::Equatorial, ReferenceFrame::Equatorial);
        assert_eq!(ReferenceFrame::Ecliptic, ReferenceFrame::Ecliptic);
        assert_eq!(ReferenceFrame::Galactic, ReferenceFrame::Galactic);
    }
}
