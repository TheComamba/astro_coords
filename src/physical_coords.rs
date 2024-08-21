//! Defines the PhysicalCartesian struct.

use std::fmt::Display;

use crate::{angle_helper::QUARTER_CIRC, reference_frame::ReferenceFrame, traits::*};

/// A wrapper around a Cartesian coordinate that is in a physical reference frame.
#[derive(Debug, Clone)]
pub struct PhysicalCoords<T>
where
    T: Mathematical + ActiveRotation<T> + AsRef<T> + Display,
{
    reference_frame: ReferenceFrame,
    mathematical_coordinates: T,
}

impl<T> PhysicalCoords<T>
where
    T: Mathematical + ActiveRotation<T> + AsRef<T> + Display,
{
    /// Create a new PhysicalCartesian object.
    ///
    /// # Example
    /// ```
    /// use astro_coords::physical_coords::PhysicalCoords;
    /// use astro_coords::reference_frame::ReferenceFrame;
    /// use astro_coords::direction::Direction;
    ///
    /// let physical = PhysicalCoords::new(Direction::X, ReferenceFrame::Equatorial);
    /// println!("{}", physical);
    /// ```
    pub fn new(inner: T, reference_frame: ReferenceFrame) -> Self {
        Self {
            reference_frame,
            mathematical_coordinates: inner,
        }
    }
}

impl<T> Physical<T> for PhysicalCoords<T>
where
    T: Mathematical + ActiveRotation<T> + AsRef<T> + Clone + Display,
{
    /// Returns the frame of reference that the mathematical coordinates are defined in.
    ///
    /// # Example
    /// ```
    /// use astro_coords::physical_coords::PhysicalCoords;
    /// use astro_coords::reference_frame::ReferenceFrame;
    /// use astro_coords::direction::Direction;
    /// use astro_coords::traits::*;
    ///
    /// let physical = PhysicalCoords::new(Direction::X, ReferenceFrame::Equatorial);
    /// assert_eq!(physical.reference_frame(), ReferenceFrame::Equatorial);
    /// ```
    fn reference_frame(&self) -> ReferenceFrame {
        self.reference_frame
    }

    /// Changes the frame of reference, transforming the mathematical coordinates.
    ///
    /// See the `in_reference_frame` method for a more detailed explanation of the transformation as well as a more detailed example.
    ///
    /// # Example
    /// ```
    /// use astro_coords::physical_coords::PhysicalCoords;
    /// use astro_coords::reference_frame::ReferenceFrame;
    /// use astro_coords::direction::Direction;
    /// use astro_coords::traits::*;
    ///
    /// let mut physical = PhysicalCoords::new(Direction::Y, ReferenceFrame::Equatorial);
    /// physical.change_reference_frame(ReferenceFrame::Ecliptic);
    /// ```
    fn change_reference_frame(&mut self, new_frame: ReferenceFrame) {
        if self.reference_frame == new_frame {
            return;
        }

        if new_frame != ReferenceFrame::Equatorial {
            if self.reference_frame != ReferenceFrame::Equatorial {
                self.change_reference_frame(ReferenceFrame::Equatorial);
            }
            let new_z = new_frame.z_axis();
            let equinox_to_q = QUARTER_CIRC - new_z.longitude;
            let plane_tilt = QUARTER_CIRC - new_z.latitude;
            let w = new_frame.prime_meridian_offset();
            self.mathematical_coordinates = self
                .mathematical_coordinates
                .rotated_z(equinox_to_q)
                .rotated_x(plane_tilt)
                .rotated_z(w);
        } else {
            let old_z = self.reference_frame.z_axis();
            let equinox_to_q = QUARTER_CIRC - old_z.longitude;
            let plane_tilt = QUARTER_CIRC - old_z.latitude;
            let w = self.reference_frame.prime_meridian_offset();
            self.mathematical_coordinates = self
                .mathematical_coordinates
                .rotated_z(-w)
                .rotated_x(-plane_tilt)
                .rotated_z(-equinox_to_q);
        }
        self.reference_frame = new_frame;
    }

    /// Returns a new instance with the mathematical coordinates transformed to the new frame of reference.
    ///
    /// The transformation is defined by a new z-axis given in right ascension and declination, and the prime meridian offset that specifies the orientation of the new x-axis.
    /// Based on Fig. 1 of [the IAU report](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf), the algorithm is as follows:
    /// - If the coordinates are not currently given in equatorial coordinates, convert them to equatorial coordinates by following the inverse steps of this algorithm.
    /// - Rotate around the old z-axis by 90°-a, where a is the right ascension of the new z-axis. The current x-axis now points to the intermediate direction Q.
    /// - Rotate around Q by 90°-d, where d is the declination of the new z-axis. The current z-axis now points to the new z-axis.
    /// - Rotate around the new z-axis by W, the prime meridian offset. The current x-axis now points to the new x-axis.
    ///
    /// # Example
    /// ```
    /// use astro_coords::physical_coords::PhysicalCoords;
    /// use astro_coords::reference_frame::ReferenceFrame;
    /// use astro_coords::spherical::Spherical;
    /// use astro_coords::traits::*;
    /// use simple_si_units::geometry::Angle;
    ///
    /// // As an example, the position of the star Sirius is expressed in various reference frames.
    /// // The values are taken from [NASA's HEASARC Object Position Finder Tool](https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/convcoord/convcoord.pl?CoordVal=Sirius&CoordType=J2000&Resolver=GRB%2FSIMBAD%2BSesame%2FNED&NoCache=on&Epoch=)
    /// let equatorial_lon = Angle::from_deg(101.287155);
    /// let equatorial_lat = Angle::from_deg(-16.716116);
    /// let equatorial = PhysicalCoords::new(Spherical::new(equatorial_lon, equatorial_lat), ReferenceFrame::Equatorial);
    /// let ecliptic_lon = Angle::from_deg(104.081665);
    /// let ecliptic_lat = Angle::from_deg(-39.605249);
    /// let ecliptic = PhysicalCoords::new(Spherical::new(ecliptic_lon, ecliptic_lat), ReferenceFrame::Ecliptic);
    /// let galactic_lon = Angle::from_deg(227.230283);
    /// let galactic_lat = Angle::from_deg(-8.890284);
    /// let galactic = PhysicalCoords::new(Spherical::new(galactic_lon, galactic_lat), ReferenceFrame::Galactic);
    ///
    /// let equatorial_from_ecliptic = ecliptic.in_reference_frame(ReferenceFrame::Equatorial);
    /// let equatorial_from_galactic = galactic.in_reference_frame(ReferenceFrame::Equatorial);
    /// let ecliptic_from_equatorial = equatorial.in_reference_frame(ReferenceFrame::Ecliptic);
    /// let ecliptic_from_galactic = galactic.in_reference_frame(ReferenceFrame::Ecliptic);
    /// let galactic_from_equatorial = equatorial.in_reference_frame(ReferenceFrame::Galactic);
    /// let galactic_from_ecliptic = ecliptic.in_reference_frame(ReferenceFrame::Galactic);
    ///
    /// let acc = Angle::from_deg(1e-5);
    /// assert!(equatorial_from_ecliptic.mathematical_coordinates().eq_within(equatorial.mathematical_coordinates(), acc));
    /// assert!(equatorial_from_galactic.mathematical_coordinates().eq_within(equatorial.mathematical_coordinates(), acc));
    /// assert!(ecliptic_from_equatorial.mathematical_coordinates().eq_within(ecliptic.mathematical_coordinates(), acc));
    /// assert!(ecliptic_from_galactic.mathematical_coordinates().eq_within(ecliptic.mathematical_coordinates(), acc));
    /// assert!(galactic_from_equatorial.mathematical_coordinates().eq_within(galactic.mathematical_coordinates(), acc));
    /// assert!(galactic_from_ecliptic.mathematical_coordinates().eq_within(galactic.mathematical_coordinates(), acc));
    /// ```
    fn in_reference_frame(&self, new_frame: ReferenceFrame) -> PhysicalCoords<T> {
        let mut new = self.clone();
        new.change_reference_frame(new_frame);
        new
    }

    /// Overwrites the frame of reference without transforming the mathematical coordinates.
    ///
    /// # Example
    /// ```
    /// use astro_coords::physical_coords::PhysicalCoords;
    /// use astro_coords::reference_frame::ReferenceFrame;
    /// use astro_coords::direction::Direction;
    /// use astro_coords::traits::*;
    ///
    /// let mut physical = PhysicalCoords::new(Direction::X, ReferenceFrame::Equatorial);
    /// physical.overwrite_reference_frame(ReferenceFrame::Ecliptic);
    /// assert_eq!(physical.reference_frame(), ReferenceFrame::Ecliptic);
    /// ```
    fn overwrite_reference_frame(&mut self, new_frame: ReferenceFrame) {
        self.reference_frame = new_frame;
    }

    /// Returns a reference to the mathematical coordinates object as it is represented in the physical frame of reference.
    ///
    /// # Example
    /// ```
    /// use astro_coords::physical_coords::PhysicalCoords;
    /// use astro_coords::reference_frame::ReferenceFrame;
    /// use astro_coords::direction::Direction;
    /// use astro_coords::traits::*;
    ///
    /// let physical = PhysicalCoords::new(Direction::X, ReferenceFrame::Equatorial);
    /// assert_eq!(physical.mathematical_coordinates(), &Direction::X);
    /// ```
    fn mathematical_coordinates(&self) -> &T {
        &self.mathematical_coordinates
    }
}

impl<T> ActiveRotation<PhysicalCoords<T>> for PhysicalCoords<T>
where
    T: Mathematical + ActiveRotation<T> + AsRef<T> + Display,
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

impl Display for PhysicalCoords<crate::direction::Direction> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} in {} coordinates",
            self.mathematical_coordinates, self.reference_frame
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::direction::Direction;

    use super::*;

    #[test]
    fn changing_from_equatorial_to_ecliptic_preserves_x_axis() {
        let equatorial = ReferenceFrame::Equatorial;
        let ecliptic = ReferenceFrame::Ecliptic;

        let mut physical = PhysicalCoords::new(Direction::X, equatorial);
        physical.change_reference_frame(ecliptic);

        assert!(physical
            .mathematical_coordinates
            .eq_within(&Direction::X, 1e-5));
    }
}
