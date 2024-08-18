//! Defines the PhysicalCartesian struct.

use crate::{reference_frame::ReferenceFrame, traits::*};

/// A wrapper around a Cartesian coordinate that is in a physical reference frame.
pub struct PhysicalCoords<T>
where
    T: Mathematical + ActiveRotation<T>,
{
    reference_frame: ReferenceFrame,
    mathematical_coordinates: T,
}

impl<T> PhysicalCoords<T>
where
    T: Mathematical + ActiveRotation<T>,
{
    /// Create a new PhysicalCartesian.
    pub fn new(inner: T, reference_frame: ReferenceFrame) -> Self {
        Self {
            reference_frame,
            mathematical_coordinates: inner,
        }
    }
}

impl<T> Physical for PhysicalCoords<T>
where
    T: Mathematical + ActiveRotation<T>,
{
    fn reference_frame(&self) -> ReferenceFrame {
        self.reference_frame
    }
}

impl<T> ActiveRotation<PhysicalCoords<T>> for PhysicalCoords<T>
where
    T: Mathematical + ActiveRotation<T>,
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
