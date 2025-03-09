//! This module contains the EarthEquatorial struct and its implementation.

use std::fmt::Display;

use uom::si::f64::{Angle, Length};

use crate::{
    angle_helper::EARTH_AXIS_TILT, cartesian::Cartesian, equatorial::Equatorial, traits::*,
};

use super::{direction::Direction, ecliptic::Ecliptic, spherical::Spherical};

/// The EarthEquatorial coordinate struct.
///
/// Earth equatorial coordinates are oriented such that the z-axis points along the Earth's axis of rotation or north pole. The remaining orientation ambiguity is resolved by any of the four following facts:
/// - At the vernal equinox, the x-axis points directly towards the sun.
/// - At the summer solstice, when the sun appears high in the sky on the northern hemisphere, its y-projection and z-projection are both positive.
/// - At the autumn equinox, the x-axis points directly away from the sun.
/// - At the winter solstice, when the sun appears high in the sky on the southern hemisphere, its y-projection and z-projection are both negative.
///
/// Note that the orientation of the coordinate system does not change with time, even though the Earth is rotating.
///
/// The coordinates in this system are expressed in right ascension and declination, a form of spherical representation.
/// - Right ascension is the angle between the vernal equinox and the point in the sky, measured along the equator.
/// - Declination is the angle between the point in the sky and the equator.
///
/// The EarthEquatorial struct provides methods to convert to and from other coordinate systems.
pub struct EarthEquatorial {
    right_ascension: Angle,
    declination: Angle,
}

impl EarthEquatorial {
    pub const fn new(right_ascension: Angle, declination: Angle) -> EarthEquatorial {
        EarthEquatorial {
            right_ascension,
            declination,
        }
    }

    pub fn to_cartesian(&self, length: Length) -> Cartesian {
        self.to_direction().to_cartesian(length)
    }

    pub fn to_direction(&self) -> Direction {
        let direction_in_equatorial =
            Spherical::new(self.right_ascension, self.declination).to_direction();
        direction_in_equatorial.rotated(-EARTH_AXIS_TILT, &Direction::X)
    }

    pub fn to_ecliptic(&self) -> Ecliptic {
        self.to_direction().to_ecliptic()
    }

    pub fn to_equatorial(&self, axis: Direction) -> Equatorial {
        self.to_direction().to_equatorial(axis)
    }

    pub fn to_spherical(&self) -> Spherical {
        self.to_direction().to_spherical()
    }

    pub fn eq_within(&self, other: &EarthEquatorial, accuracy: Angle) -> bool {
        let mut self_spherical = Spherical::new(self.right_ascension, self.declination);
        self_spherical.normalize();
        let mut other_spherical = Spherical::new(other.right_ascension, other.declination);
        other_spherical.normalize();
        self_spherical.eq_within(&other_spherical, accuracy)
    }
}

impl Display for EarthEquatorial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "RA: {}, Dec: {}", self.right_ascension, self.declination)
    }
}

#[cfg(test)]
pub(super) mod tests {
    use uom::si::f64::Angle;

    use crate::{angle_helper::*, direction::Direction, ecliptic::Ecliptic, spherical::Spherical};

    use super::EarthEquatorial;

    const TEST_ACCURACY: f64 = 1e-5;
    const ANGLE_TEST_ACCURACY: Angle = Angle { rad: TEST_ACCURACY };

    pub(super) const EARTH_NORTH_POLE_IN_ECLIPTIC_COORDINATES: Ecliptic =
        Ecliptic::new(Spherical::new(
            QUARTER_CIRC,
            Angle {
                rad: QUARTER_CIRC.get::<radian>() - EARTH_AXIS_TILT.get::<radian>(),
            },
        ));

    #[test]
    fn the_sun_is_high_in_northern_hemisphere_during_summer_solstice() {
        let direction = Direction::Y;
        let equatorial = direction.to_earth_equatorial();
        assert!(equatorial.declination.get::<radian>() > 0.);
    }

    /*
     * https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=0&dec=90&pa=0.0&out_csys=Ecliptic&out_equinox=J2000.0
     */
    #[test]
    fn ra_zero_dec_ninty_is_north_pole() {
        let equatorial = EarthEquatorial::new(Angle::new::<degree>(0.), Angle::new::<degree>(90.));
        let expected = EARTH_NORTH_POLE_IN_ECLIPTIC_COORDINATES;
        let actual = equatorial.to_ecliptic();
        println!("expected: {},\n  actual: {}", expected, actual);
        assert!(actual.eq_within(
            &EARTH_NORTH_POLE_IN_ECLIPTIC_COORDINATES,
            ANGLE_TEST_ACCURACY
        ));
    }

    /*
     * https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=0&dec=0&pa=0.0&out_csys=Ecliptic&out_equinox=J2000.0
     */
    #[test]
    fn ra_zero_dec_zer_is_x_axis() {
        let equatorial = EarthEquatorial::new(ANGLE_ZERO, ANGLE_ZERO);
        let expected = Ecliptic::X_DIRECTION;
        let actual = equatorial.to_ecliptic();
        println!("expected: {},\n  actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    /*
     * https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=6&dec=0&pa=0.0&out_csys=Ecliptic&out_equinox=J2000.0
     */
    #[test]
    fn ra_ninty_dec_zero_is_equator_zenith() {
        let equatorial = EarthEquatorial::new(Angle::new::<degree>(90.), ANGLE_ZERO);
        let expected = Ecliptic::new(Spherical::new(Angle::new::<degree>(90.), -EARTH_AXIS_TILT));
        let actual = equatorial.to_ecliptic();
        println!("expected: {},\n  actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    /*
     * https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=12&dec=0&pa=0.0&out_csys=Ecliptic&out_equinox=J2000.0
     */
    #[test]
    fn ra_oneeighty_dec_zero_is_minus_x_axis() {
        let equatorial = EarthEquatorial::new(Angle::new::<degree>(180.), ANGLE_ZERO);
        let expected = -Ecliptic::X_DIRECTION;
        let actual = equatorial.to_ecliptic();
        println!("expected: {},\n  actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    /*
     * https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=18&dec=0&pa=0.0&out_csys=Ecliptic&out_equinox=J2000.0
     */
    #[test]
    fn ra_twoseventy_dec_zero_is_equator_midnight() {
        let equatorial = EarthEquatorial::new(Angle::new::<degree>(270.), ANGLE_ZERO);
        let expected = Ecliptic::new(Spherical::new(Angle::new::<degree>(270.), EARTH_AXIS_TILT));
        let actual = equatorial.to_ecliptic();
        println!("expected: {},\n  actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    /*
     * https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=6&dec=23%2027&pa=0.0&out_csys=Ecliptic&out_equinox=J2000.0
     */
    #[test]
    fn ra_ninty_dec_tilt_is_y_axis() {
        let equatorial = EarthEquatorial::new(Angle::new::<degree>(90.), EARTH_AXIS_TILT);
        let expected = Ecliptic::Y_DIRECTION;
        let actual = equatorial.to_ecliptic();
        println!("expected: {},\n  actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    /*
     * https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=18&dec=66%2033&pa=0.0&out_csys=Ecliptic&out_equinox=J2000.0
     */
    #[test]
    fn ra_twoseventy_dec_ninty_minus_tilt_is_z_axis() {
        let equatorial = EarthEquatorial::new(
            Angle::new::<degree>(270.),
            Angle::new::<degree>(90.) - EARTH_AXIS_TILT,
        );
        let expected = Ecliptic::Z_DIRECTION;
        let actual = equatorial.to_ecliptic();
        println!("expected: {},\n  actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    /*
     * Calculated using https://frostydrew.org/utilities.dc/convert/tool-eq_coordinates/
     */
    #[test]
    fn specific_testcase() {
        let equatorial =
            EarthEquatorial::new(Angle::new::<degree>(234.), Angle::new::<degree>(56.));
        let expected = Spherical::new(
            Angle::new::<degree>(194.547656),
            Angle::new::<degree>(70.149178),
        )
        .to_ecliptic();
        let actual = equatorial.to_ecliptic();
        println!("expected: {}, actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    #[test]
    fn roundtrip() {
        const STEPS: usize = 100;
        for i in 0..STEPS {
            for j in 0..STEPS {
                let ra = Angle::new::<degree>(360. * i as f64 / STEPS as f64);
                let dec = Angle::new::<degree>(180. * j as f64 / STEPS as f64 - 90.);
                let equatorial = EarthEquatorial::new(ra, dec);
                let ecliptic_dir = equatorial.to_direction();
                let actual = ecliptic_dir.to_earth_equatorial();
                println!("expected: {}, actual: {}", equatorial, actual);
                assert!(actual.eq_within(&equatorial, ANGLE_TEST_ACCURACY));
            }
        }
    }
}
