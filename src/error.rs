//! Error types for the astro_coordinates crate

use std::fmt::Display;

#[derive(Debug, Clone, Copy)]
pub enum AstroCoordsError {
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
