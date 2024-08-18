//! Defines the PhysicalCartesian struct.

use crate::{reference_frame::ReferenceFrame, traits::*};

/// A wrapper around a Cartesian coordinate that is in a physical reference frame.
pub struct PhysicalCartesian<T>
where
    T: Mathematical + ActiveRotation<T>,
{
    reference_frame: ReferenceFrame,
    mathematical_coordinates: T,
}

impl<T> PhysicalCartesian<T>
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

impl<T> Physical for PhysicalCartesian<T>
where
    T: Mathematical + ActiveRotation<T>,
{
    fn reference_frame(&self) -> ReferenceFrame {
        self.reference_frame
    }
}
