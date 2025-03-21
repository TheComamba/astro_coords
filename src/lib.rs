#![warn(clippy::unwrap_used)]
// #![warn(missing_docs)]
#![doc = include_str!("../README.md")]
#[macro_use]
extern crate uom;

mod angle_helper;
pub mod cartesian;
pub mod direction;
pub mod earth_equatorial;
pub mod ecliptic;
pub mod equatorial;
pub mod error;
pub mod physical_coords;
pub mod ra_and_dec;
pub mod reference_frame;
pub mod spherical;
pub mod traits;
pub mod transformations;
pub mod units;

pub(crate) const NORMALIZATION_THRESHOLD: f64 = 1e-12; // When subtracting two f64 numbers around unity, the uncertainty is around 1e-17.
