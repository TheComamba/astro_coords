//! Error types for the astro_coords crate

use std::fmt::Display;

/// Error type for the AstroCoords struct
#[derive(Debug, Clone, Copy)]
pub enum AstroCoordsError {
    /// Error for trying to normalize a zero vector
    NormalizingZeroVector,
}

impl Display for AstroCoordsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AstroCoordsError::NormalizingZeroVector => {
                write!(f, "Cannot normalize a zero vector")
            }
        }
    }
}
