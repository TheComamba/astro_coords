//! Defines the PhysicalCartesian struct.

use crate::{reference_frame::ReferenceFrame, traits::*};

/// A wrapper around a Cartesian coordinate that is in a physical reference frame.
#[derive(Debug, Clone)]
pub struct PhysicalCoords<T>
where
    T: Mathematical + ActiveRotation<T> + AsRef<T>,
{
    reference_frame: ReferenceFrame,
    mathematical_coordinates: T,
}

impl<T> PhysicalCoords<T>
where
    T: Mathematical + ActiveRotation<T> + AsRef<T>,
{
    /// Create a new PhysicalCartesian object.
    pub fn new(inner: T, reference_frame: ReferenceFrame) -> Self {
        Self {
            reference_frame,
            mathematical_coordinates: inner,
        }
    }
}

impl<T> Physical for PhysicalCoords<T>
where
    T: Mathematical + ActiveRotation<T> + AsRef<T>,
{
    /// Returns the frame of reference that the mathematical coordinates are defined in.
    fn reference_frame(&self) -> ReferenceFrame {
        self.reference_frame
    }

    /// Changes the frame of reference, transforming the mathematical coordinates.
    ///
    /// The transformation is defined by a new z-axis, and an angle-parameter W that specifies the orientation of the new x-axis.
    /// Based on Fig. 1 of [the IAU report](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf), the algorithm is as follows:
    /// - Express the new z-axis in the old reference frame.
    /// - Convert that new-z-axis to spherical coordinates to determine the longitude (alpha) and latitude (delta) of the new z-axis.
    /// - Rotate around the old z-axis by 90°-a. The current x-axis now points to the intermediate direction Q.
    /// - Rotate around Q by 90°-d. The current z-axis now points to the new z-axis.
    /// - Rotate around the new z-axis by W. The current x-axis now points to the new x-axis.
    ///
    /// # Example
    /// TODO
    fn change_reference_frame(&mut self, new_frame: ReferenceFrame) {
        if self.reference_frame == new_frame {
            return;
        }
        todo!();
        self.reference_frame = new_frame;
    }

    /// Overwrites the frame of reference without transforming the mathematical coordinates.
    fn overwrite_reference_frame(&mut self, new_frame: ReferenceFrame) {
        self.reference_frame = new_frame;
    }

    fn mathematical_coordinates(&self) -> &dyn Mathematical {
        self.mathematical_coordinates.as_ref()
    }
}

impl<T> ActiveRotation<PhysicalCoords<T>> for PhysicalCoords<T>
where
    T: Mathematical + ActiveRotation<T> + AsRef<T>,
{
    fn rotated(
        &self,
        angle: simple_si_units::geometry::Angle<f64>,
        axis: &crate::direction::Direction,
    ) -> PhysicalCoords<T> {
        Self {
            reference_frame: self.reference_frame,
            mathematical_coordinates: self.mathematical_coordinates.rotated(angle, axis),
        }
    }

    fn rotated_x(&self, angle: simple_si_units::geometry::Angle<f64>) -> PhysicalCoords<T> {
        Self {
            reference_frame: self.reference_frame,
            mathematical_coordinates: self.mathematical_coordinates.rotated_x(angle),
        }
    }

    fn rotated_y(&self, angle: simple_si_units::geometry::Angle<f64>) -> PhysicalCoords<T> {
        Self {
            reference_frame: self.reference_frame,
            mathematical_coordinates: self.mathematical_coordinates.rotated_y(angle),
        }
    }

    fn rotated_z(&self, angle: simple_si_units::geometry::Angle<f64>) -> PhysicalCoords<T> {
        Self {
            reference_frame: self.reference_frame,
            mathematical_coordinates: self.mathematical_coordinates.rotated_z(angle),
        }
    }

    fn active_rotation_to_new_z_axis(&self, new_z: &PhysicalCoords<T>) -> PhysicalCoords<T> {
        Self {
            reference_frame: self.reference_frame,
            mathematical_coordinates: self
                .mathematical_coordinates
                .active_rotation_to_new_z_axis(&new_z.mathematical_coordinates),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::direction::Direction;

    use super::*;

    #[test]
    fn changing_from_equatorial_to_ecliptic_preserves_x_axis() {
        let equatorial = ReferenceFrame::Equatorial;
        let ecliptic = ReferenceFrame::Ecliptic;

        let mut physical = PhysicalCoords::new(Direction::X, equatorial);
        physical.change_reference_frame(ecliptic);

        assert!(physical
            .mathematical_coordinates
            .eq_within(&Direction::X, 1e-5));
    }
}
