//! Defines the PhysicalCartesian struct.

use crate::{cartesian::Cartesian, reference_frame::ReferenceFrame};

/// A wrapper around a Cartesian coordinate that is in a physical reference frame.
pub struct PhysicalCartesian {
    reference_frame: ReferenceFrame,
    inner: Cartesian,
}

impl PhysicalCartesian {
    /// Create a new PhysicalCartesian.
    pub fn new(inner: Cartesian, reference_frame: ReferenceFrame) -> Self {
        Self {
            reference_frame,
            inner,
        }
    }

    /// Get the reference frame.
    pub fn reference_frame(&self) -> ReferenceFrame {
        self.reference_frame
    }
}
