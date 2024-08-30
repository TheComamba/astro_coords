/// A marker trait for types that represent mathematical coordinates.
///
/// This trait is used to ensure that the coordinates do not have any pre-defined physical meaning. More specifically, coordinates implementing this trait lack a frame of reference.
pub trait Mathematical {}
