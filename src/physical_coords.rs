//! Defines the PhysicalCartesian struct.

use std::fmt::Display;

use uom::si::f64::{Angle, Time};

use crate::{
    angle_helper::*,
    reference_frame::{CelestialBody, ReferenceFrame},
    traits::*,
};

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

    fn change_from_ecliptic_to_equatorial(&mut self) {
        self.mathematical_coordinates = self.mathematical_coordinates.rotated_x(EARTH_AXIS_TILT());
    }

    fn change_from_galactic_to_equatorial(&mut self) {
        self.mathematical_coordinates = self
            .mathematical_coordinates
            .rotated_z(HALF_CIRC() - GALACTIC_LONGITUDE_OF_NORTH_CELESTIAL_POLE())
            .rotated_y(QUARTER_CIRC() - DEC_OF_GALACTIC_NORTH())
            .rotated_z(RA_OF_GALACTIC_NORTH());
    }

    fn change_from_cartographic_to_equatorial(
        &mut self,
        celestial_body: CelestialBody,
        time_since_epoch: Time,
    ) {
        let old_z = celestial_body.z_axis();
        let equinox_to_q = QUARTER_CIRC() - old_z.longitude;
        let plane_tilt = QUARTER_CIRC() - old_z.latitude;
        let w = celestial_body.prime_meridian_offset(time_since_epoch);
        self.mathematical_coordinates = self
            .mathematical_coordinates
            .rotated_z(-w)
            .rotated_x(-plane_tilt)
            .rotated_z(-equinox_to_q);
    }

    fn change_to_equatorial(&mut self) {
        match self.reference_frame {
            ReferenceFrame::Equatorial => {}
            ReferenceFrame::Ecliptic => self.change_from_ecliptic_to_equatorial(),
            ReferenceFrame::Galactic => self.change_from_galactic_to_equatorial(),
            ReferenceFrame::Cartographic(body, time) => {
                self.change_from_cartographic_to_equatorial(body, time)
            }
        }
    }

    fn change_from_equatorial_to_ecliptic(&mut self) {
        self.mathematical_coordinates = self.mathematical_coordinates.rotated_x(-EARTH_AXIS_TILT());
    }

    fn change_from_equatorial_to_galactic(&mut self) {
        self.mathematical_coordinates = self
            .mathematical_coordinates
            .rotated_z(-RA_OF_GALACTIC_NORTH())
            .rotated_y(DEC_OF_GALACTIC_NORTH() - QUARTER_CIRC())
            .rotated_z(GALACTIC_LONGITUDE_OF_NORTH_CELESTIAL_POLE() - HALF_CIRC());
    }

    fn change_from_equatorial_to_cartographic(
        &mut self,
        celestial_body: CelestialBody,
        time_since_epoch: Time,
    ) {
        let new_z = celestial_body.z_axis();
        let equinox_to_q = QUARTER_CIRC() - new_z.longitude;
        let plane_tilt = QUARTER_CIRC() - new_z.latitude;
        let w = celestial_body.prime_meridian_offset(time_since_epoch);
        self.mathematical_coordinates = self
            .mathematical_coordinates
            .rotated_z(equinox_to_q)
            .rotated_x(plane_tilt)
            .rotated_z(w);
    }

    fn change_from_equatorial(&mut self, new_frame: ReferenceFrame) {
        match new_frame {
            ReferenceFrame::Equatorial => {}
            ReferenceFrame::Ecliptic => self.change_from_equatorial_to_ecliptic(),
            ReferenceFrame::Galactic => self.change_from_equatorial_to_galactic(),
            ReferenceFrame::Cartographic(body, time) => {
                self.change_from_equatorial_to_cartographic(body, time)
            }
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
                self.change_to_equatorial();
            }
            self.change_from_equatorial(new_frame);
        } else {
            self.change_to_equatorial();
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
    /// use uom::si::f64::Angle;
    ///
    /// // As an example, the position of the star Sirius is expressed in various reference frames.
    /// // The values are taken from [NASA's HEASARC Object Position Finder Tool](https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/convcoord/convcoord.pl?CoordVal=Sirius&CoordType=J2000&Resolver=GRB%2FSIMBAD%2BSesame%2FNED&NoCache=on&Epoch=)
    /// let equatorial_lon = Angle::new::<degree>(101.287155);
    /// let equatorial_lat = Angle::new::<degree>(-16.716116);
    /// let equatorial = PhysicalCoords::new(Spherical::new(equatorial_lon, equatorial_lat), ReferenceFrame::Equatorial);
    /// let ecliptic_lon = Angle::new::<degree>(104.081665);
    /// let ecliptic_lat = Angle::new::<degree>(-39.605249);
    /// let ecliptic = PhysicalCoords::new(Spherical::new(ecliptic_lon, ecliptic_lat), ReferenceFrame::Ecliptic);
    /// let galactic_lon = Angle::new::<degree>(227.230283);
    /// let galactic_lat = Angle::new::<degree>(-8.890284);
    /// let galactic = PhysicalCoords::new(Spherical::new(galactic_lon, galactic_lat), ReferenceFrame::Galactic);
    ///
    /// let equatorial_from_ecliptic = ecliptic.in_reference_frame(ReferenceFrame::Equatorial);
    /// let equatorial_from_galactic = galactic.in_reference_frame(ReferenceFrame::Equatorial);
    /// let ecliptic_from_equatorial = equatorial.in_reference_frame(ReferenceFrame::Ecliptic);
    /// let ecliptic_from_galactic = galactic.in_reference_frame(ReferenceFrame::Ecliptic);
    /// let galactic_from_equatorial = equatorial.in_reference_frame(ReferenceFrame::Galactic);
    /// let galactic_from_ecliptic = ecliptic.in_reference_frame(ReferenceFrame::Galactic);
    ///
    /// let acc = Angle::new::<degree>(1e-1); // Rotations can include some significant loss of accuracy.
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
    fn rotated(&self, angle: Angle, axis: &crate::direction::Direction) -> PhysicalCoords<T> {
        Self {
            reference_frame: self.reference_frame,
            mathematical_coordinates: self.mathematical_coordinates.rotated(angle, axis),
        }
    }

    fn rotated_x(&self, angle: Angle) -> PhysicalCoords<T> {
        Self {
            reference_frame: self.reference_frame,
            mathematical_coordinates: self.mathematical_coordinates.rotated_x(angle),
        }
    }

    fn rotated_y(&self, angle: Angle) -> PhysicalCoords<T> {
        Self {
            reference_frame: self.reference_frame,
            mathematical_coordinates: self.mathematical_coordinates.rotated_y(angle),
        }
    }

    fn rotated_z(&self, angle: Angle) -> PhysicalCoords<T> {
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
    use std::f64::consts::PI;

    use uom::si::angle::{degree, radian};

    use crate::{
        angle_helper::{normalized_angle, EARTH_AXIS_TILT},
        direction::Direction,
        ra_and_dec::{Declination, RightAscension, Sgn},
        spherical::Spherical,
    };

    use super::*;

    #[test]
    fn changing_from_equatorial_to_ecliptic_preserves_x_axis() {
        let equatorial = ReferenceFrame::Equatorial;
        let ecliptic = ReferenceFrame::Ecliptic;

        let physical = PhysicalCoords::new(Direction::X, equatorial).in_reference_frame(ecliptic);

        assert!(physical
            .mathematical_coordinates
            .eq_within(&Direction::X, 1e-5));
    }

    #[test]
    fn the_equations_for_changing_from_equatorial_to_ecliptic_reference_hold() {
        // https://aas.aanda.org/articles/aas/full/1998/01/ds1449/node3.html

        let e = EARTH_AXIS_TILT().get::<radian>();

        let angles = vec![0., PI, 0.3, 1.4, -1.4, 3.5, 7.];
        for ra in angles.clone() {
            for dec in angles.clone() {
                let equatorial = PhysicalCoords::new(
                    Spherical::new(Angle::new::<radian>(ra), Angle::new::<radian>(dec)),
                    ReferenceFrame::Equatorial,
                );
                let ecliptic = equatorial.in_reference_frame(ReferenceFrame::Ecliptic);
                let ecliptic_coords = ecliptic.mathematical_coordinates();
                let l = ecliptic_coords.longitude.get::<radian>();
                let b = ecliptic_coords.latitude.get::<radian>();

                let lefthand = b.sin();
                let righthand = dec.sin() * e.cos() - dec.cos() * e.sin() * ra.sin();
                assert!(
                    (lefthand - righthand).abs() < 1e-5,
                    "ra: {}, dec: {}",
                    ra,
                    dec
                );

                let lefthand = ra.cos() * dec.cos();
                let righthand = l.cos() * b.cos();
                assert!(
                    (lefthand - righthand).abs() < 1e-5,
                    "ra: {}, dec: {}",
                    ra,
                    dec
                );

                let lefthand = l.sin() * b.cos();
                let righthand = dec.sin() * e.sin() + dec.cos() * e.cos() * ra.sin();
                assert!(
                    (lefthand - righthand).abs() < 1e-5,
                    "ra: {}, dec: {}",
                    ra,
                    dec
                );
            }
        }
    }

    #[test]
    fn the_equations_for_changing_from_equatorial_to_galactic_reference_hold() {
        // https://en.wikipedia.org/wiki/Galactic_coordinate_system#Conversion_between_equatorial_and_galactic_coordinates
        let angp = RightAscension::new(12, 51, 24.).to_angle().get::<radian>() as f64;
        let dngp = Angle::new::<degree>(27.13).get::<radian>() as f64;
        let lncp = Angle::new::<degree>(122.93314).get::<radian>() as f64;

        let angles = vec![0., PI, 0.3, 1.4, -1.4, 3.5, 7.];
        for ra in angles.clone() {
            for dec in angles.clone() {
                let equatorial = PhysicalCoords::new(
                    Spherical::new(Angle::new::<radian>(ra), Angle::new::<radian>(dec)),
                    ReferenceFrame::Equatorial,
                );
                let galactic = equatorial.in_reference_frame(ReferenceFrame::Galactic);
                let galactic_coords = galactic.mathematical_coordinates();
                let l = galactic_coords.longitude.get::<radian>();
                let b = galactic_coords.latitude.get::<radian>();

                let lefthand = b.sin();
                let righthand = dngp.sin() * dec.sin() + dngp.cos() * dec.cos() * (ra - angp).cos();
                assert!(
                    (lefthand - righthand).abs() < 1e-5,
                    "ra: {}, dec: {}",
                    ra,
                    dec
                );

                let lefthand = b.cos() * (lncp - l).sin();
                let righthand = dec.cos() * (ra - angp).sin();
                assert!(
                    (lefthand - righthand).abs() < 1e-5,
                    "ra: {}, dec: {}",
                    ra,
                    dec
                );

                let lefthand = b.cos() * (lncp - l).cos();
                let righthand = dngp.cos() * dec.sin() - dngp.sin() * dec.cos() * (ra - angp).cos();
                assert!(
                    (lefthand - righthand).abs() < 1e-5,
                    "ra: {}, dec: {}",
                    ra,
                    dec
                );
            }
        }
    }

    #[test]
    fn specific_test_case() {
        // https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=11&dec=22%2033&pa=0.0&out_csys=Equatorial&out_equinox=J2000.0
        let ra = RightAscension::new(11, 0, 0.).to_angle();
        let dec = Declination::new(Sgn::Pos, 22, 33, 0.).to_angle();
        let equatorial = PhysicalCoords::new(Spherical::new(ra, dec), ReferenceFrame::Equatorial);

        let ecliptic = equatorial.in_reference_frame(ReferenceFrame::Ecliptic);
        let ecliptic_coords = ecliptic.mathematical_coordinates();
        let l = ecliptic_coords.longitude;
        let b = ecliptic_coords.latitude;
        let expected_l = Angle::new::<degree>(157.371833);
        let expected_b = Angle::new::<degree>(14.878115);
        assert!((l - expected_l).get::<radian>().abs() < 1e-5);
        assert!((b - expected_b).get::<radian>().abs() < 1e-5);

        let galactic = equatorial.in_reference_frame(ReferenceFrame::Galactic);
        let galactic_coords = galactic.mathematical_coordinates();
        let l = galactic_coords.longitude;
        let b = galactic_coords.latitude;
        let expected_l = normalized_angle(Angle::new::<degree>(217.042106));
        let expected_b = Angle::new::<degree>(64.361644);
        // Galactic coordinates seem to not be defined as accurately as ecliptic coordinates
        assert!((l - expected_l).get::<radian>().abs() < 2e-4);
        assert!((b - expected_b).get::<radian>().abs() < 2e-4);
    }
}
