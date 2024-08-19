use crate::reference_frame::ReferenceFrame;

pub trait Physical {
    /// Returns the frame of reference that the mathematical coordinates are defined in.
    fn reference_frame(&self) -> ReferenceFrame;

    /// Changes the frame of reference, transforming the mathematical coordinates.
    fn change_reference_frame(&mut self, new_frame: ReferenceFrame);

    /// Overwrites the frame of reference without transforming the mathematical coordinates.
    fn overwrite_reference_frame(&mut self, new_frame: ReferenceFrame);
}
