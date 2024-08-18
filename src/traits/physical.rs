use crate::reference_frame::ReferenceFrame;

pub trait Physical {
    /// Returns the frame of reference that the mathematical coordinates are defined in.
    fn reference_frame(&self) -> ReferenceFrame;
}
