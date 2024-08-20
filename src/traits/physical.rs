use crate::reference_frame::ReferenceFrame;

use super::Mathematical;

pub trait Physical<T: Mathematical> {
    /// Returns the frame of reference that the mathematical coordinates are defined in.
    fn reference_frame(&self) -> ReferenceFrame;

    /// Returns a reference to the mathematical coordinates.
    fn mathematical_coordinates(&self) -> &T;

    /// Changes the frame of reference, transforming the mathematical coordinates.
    fn change_reference_frame(&mut self, new_frame: ReferenceFrame);

    /// Overwrites the frame of reference without transforming the mathematical coordinates.
    fn overwrite_reference_frame(&mut self, new_frame: ReferenceFrame);
}
