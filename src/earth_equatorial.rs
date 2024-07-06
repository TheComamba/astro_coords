use std::fmt::Display;

use super::{direction::Direction, ecliptic::EclipticCoordinates, spherical::SphericalCoordinates};
use crate::{astro_display::AstroDisplay, real_data::planets::EARTH};
use simple_si_units::{base::Distance, geometry::Angle};

pub struct EarthEquatorialCoordinates {
    right_ascension: Angle<f64>,
    declination: Angle<f64>,
}

impl EarthEquatorialCoordinates {
    pub const fn new(
        right_ascension: Angle<f64>,
        declination: Angle<f64>,
    ) -> EarthEquatorialCoordinates {
        EarthEquatorialCoordinates {
            right_ascension,
            declination,
        }
    }

    pub(crate) fn to_direction(&self) -> Direction {
        let direction_in_equatorial =
            SphericalCoordinates::new(self.right_ascension, self.declination).to_direction();
        direction_in_equatorial.rotated(-EARTH.axis_tilt, &Direction::X)
    }

    pub fn to_ecliptic(&self) -> EclipticCoordinates {
        let dir = self.to_direction();
        let vec = dir.to_cartesian(Distance::from_meters(1.));
        vec.to_ecliptic()
    }

    pub fn to_spherical(&self) -> SphericalCoordinates {
        SphericalCoordinates::new(self.right_ascension, self.declination)
    }

    #[cfg(test)]
    pub(crate) fn eq_within(
        &self,
        other: &EarthEquatorialCoordinates,
        accuracy: Angle<f64>,
    ) -> bool {
        let mut self_spherical = SphericalCoordinates::new(self.right_ascension, self.declination);
        self_spherical.normalize();
        let mut other_spherical =
            SphericalCoordinates::new(other.right_ascension, other.declination);
        other_spherical.normalize();
        self_spherical.eq_within(&other_spherical, accuracy)
    }
}

impl AstroDisplay for EarthEquatorialCoordinates {
    fn astro_display(&self) -> String {
        format!("(ra: {}, dec: {})", self.right_ascension, self.declination)
    }
}

impl Display for EarthEquatorialCoordinates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.astro_display())
    }
}

#[cfg(test)]
pub(super) const EARTH_NORTH_POLE_IN_ECLIPTIC_COORDINATES: EclipticCoordinates =
    EclipticCoordinates::new(SphericalCoordinates::new(
        crate::units::angle::QUARTER_CIRC,
        Angle {
            rad: crate::units::angle::QUARTER_CIRC.rad - EARTH.axis_tilt.rad,
        },
    ));

#[cfg(test)]
mod tests {
    use super::EarthEquatorialCoordinates;
    use crate::{
        coordinates::{
            earth_equatorial::EARTH_NORTH_POLE_IN_ECLIPTIC_COORDINATES,
            ecliptic::EclipticCoordinates, spherical::SphericalCoordinates,
        },
        real_data::planets::*,
        units::{
            angle::{angle_eq_within, ANGLE_ZERO},
            tests::ANGLE_TEST_ACCURACY,
        },
    };
    use simple_si_units::geometry::Angle;

    const TILT_ACCURACY: Angle<f64> = Angle { rad: 2e-3 };

    /*
     * https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=0&dec=90&pa=0.0&out_csys=Ecliptic&out_equinox=J2000.0
     */
    #[test]
    fn ra_zero_dec_ninty_is_north_pole() {
        let equatorial =
            EarthEquatorialCoordinates::new(Angle::from_degrees(0.), Angle::from_degrees(90.));
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
        let equatorial = EarthEquatorialCoordinates::new(ANGLE_ZERO, ANGLE_ZERO);
        let expected = EclipticCoordinates::X_DIRECTION;
        let actual = equatorial.to_ecliptic();
        println!("expected: {},\n  actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    /*
     * https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=6&dec=0&pa=0.0&out_csys=Ecliptic&out_equinox=J2000.0
     */
    #[test]
    fn ra_ninty_dec_zero_is_equator_zenith() {
        let equatorial = EarthEquatorialCoordinates::new(Angle::from_degrees(90.), ANGLE_ZERO);
        let expected = EclipticCoordinates::new(SphericalCoordinates::new(
            Angle::from_degrees(90.),
            -EARTH.axis_tilt,
        ));
        let actual = equatorial.to_ecliptic();
        println!("expected: {},\n  actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    /*
     * https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=12&dec=0&pa=0.0&out_csys=Ecliptic&out_equinox=J2000.0
     */
    #[test]
    fn ra_oneeighty_dec_zero_is_minus_x_axis() {
        let equatorial = EarthEquatorialCoordinates::new(Angle::from_degrees(180.), ANGLE_ZERO);
        let expected = -EclipticCoordinates::X_DIRECTION;
        let actual = equatorial.to_ecliptic();
        println!("expected: {},\n  actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    /*
     * https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=18&dec=0&pa=0.0&out_csys=Ecliptic&out_equinox=J2000.0
     */
    #[test]
    fn ra_twoseventy_dec_zero_is_equator_midnight() {
        let equatorial = EarthEquatorialCoordinates::new(Angle::from_degrees(270.), ANGLE_ZERO);
        let expected = EclipticCoordinates::new(SphericalCoordinates::new(
            Angle::from_degrees(270.),
            EARTH.axis_tilt,
        ));
        let actual = equatorial.to_ecliptic();
        println!("expected: {},\n  actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    /*
     * https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=6&dec=23%2027&pa=0.0&out_csys=Ecliptic&out_equinox=J2000.0
     */
    #[test]
    fn ra_ninty_dec_tilt_is_y_axis() {
        let equatorial = EarthEquatorialCoordinates::new(Angle::from_degrees(90.), EARTH.axis_tilt);
        let expected = EclipticCoordinates::Y_DIRECTION;
        let actual = equatorial.to_ecliptic();
        println!("expected: {},\n  actual: {}", expected, actual);
        assert!(actual.eq_within(&expected, ANGLE_TEST_ACCURACY));
    }

    /*
     * https://ned.ipac.caltech.edu/coordinate_calculator?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&ra=18&dec=66%2033&pa=0.0&out_csys=Ecliptic&out_equinox=J2000.0
     */
    #[test]
    fn ra_twoseventy_dec_ninty_minus_tilt_is_z_axis() {
        let equatorial = EarthEquatorialCoordinates::new(
            Angle::from_degrees(270.),
            Angle::from_degrees(90.) - EARTH.axis_tilt,
        );
        let expected = EclipticCoordinates::Z_DIRECTION;
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
            EarthEquatorialCoordinates::new(Angle::from_degrees(234.), Angle::from_degrees(56.));
        let expected = SphericalCoordinates::new(
            Angle::from_degrees(194.547656),
            Angle::from_degrees(70.149178),
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
                let ra = Angle::from_degrees(360. * i as f64 / STEPS as f64);
                let dec = Angle::from_degrees(180. * j as f64 / STEPS as f64 - 90.);
                let equatorial = EarthEquatorialCoordinates::new(ra, dec);
                let ecliptic_dir = equatorial.to_direction();
                let actual = ecliptic_dir.to_earth_equatorial();
                println!("expected: {}, actual: {}", equatorial, actual);
                assert!(actual.eq_within(&equatorial, ANGLE_TEST_ACCURACY));
            }
        }
    }

    #[test]
    fn axis_tilt_of_mercury() {
        let orbit_normal = MERCURY.orbit.normal();
        let north = MERCURY.rotation_axis.to_direction();
        println!("orbit_normal: {}", orbit_normal);
        println!("north: {}", north);
        let expected = MERCURY.axis_tilt;
        let actual = orbit_normal.angle_to(&north);
        println!("expected: {}, actual: {}", expected, actual);
        assert!(angle_eq_within(actual, expected, TILT_ACCURACY));
    }

    #[test]
    fn axis_tilt_of_venus() {
        let orbit_normal = VENUS.orbit.normal();
        let north = VENUS.rotation_axis.to_direction();
        println!("orbit_normal: {}", orbit_normal);
        println!("north: {}", north);
        let expected = VENUS.axis_tilt;
        let actual = orbit_normal.angle_to(&north);
        println!("expected: {}, actual: {}", expected, actual);
        assert!(angle_eq_within(actual, expected, TILT_ACCURACY));
    }

    #[test]
    fn axis_tilt_of_earth() {
        let orbit_normal = EARTH.orbit.normal();
        let north = EARTH.rotation_axis.to_direction();
        println!("orbit_normal: {}", orbit_normal);
        println!("north: {}", north);
        let expected = EARTH.axis_tilt;
        let actual = orbit_normal.angle_to(&north);
        println!("expected: {}, actual: {}", expected, actual);
        assert!(angle_eq_within(actual, expected, TILT_ACCURACY));
    }

    #[test]
    fn axis_tilt_of_mars() {
        let orbit_normal = MARS.orbit.normal();
        let north = MARS.rotation_axis.to_direction();
        println!("orbit_normal: {}", orbit_normal);
        println!("north: {}", north);
        let expected = MARS.axis_tilt;
        let actual = orbit_normal.angle_to(&north);
        println!("expected: {}, actual: {}", expected, actual);
        assert!(angle_eq_within(actual, expected, TILT_ACCURACY));
    }

    #[test]
    fn axis_tilt_of_ceres() {
        let orbit_normal = CERES.orbit.normal();
        let north = CERES.rotation_axis.to_direction();
        println!("orbit_normal: {}", orbit_normal);
        println!("north: {}", north);
        let expected = CERES.axis_tilt;
        let actual = orbit_normal.angle_to(&north);
        println!("expected: {}, actual: {}", expected, actual);
        assert!(angle_eq_within(actual, expected, TILT_ACCURACY));
    }

    #[test]
    fn axis_tilt_of_jupiter() {
        let orbit_normal = JUPITER.orbit.normal();
        let north = JUPITER.rotation_axis.to_direction();
        println!("orbit_normal: {}", orbit_normal);
        println!("north: {}", north);
        let expected = JUPITER.axis_tilt;
        let actual = orbit_normal.angle_to(&north);
        println!("expected: {}, actual: {}", expected, actual);
        assert!(angle_eq_within(actual, expected, TILT_ACCURACY));
    }

    #[test]
    fn axis_tilt_of_saturn() {
        let orbit_normal = SATURN.orbit.normal();
        let north = SATURN.rotation_axis.to_direction();
        println!("orbit_normal: {}", orbit_normal);
        println!("north: {}", north);
        let expected = SATURN.axis_tilt;
        let actual = orbit_normal.angle_to(&north);
        println!("expected: {}, actual: {}", expected, actual);
        assert!(angle_eq_within(actual, expected, TILT_ACCURACY));
    }

    #[test]
    fn axis_tilt_of_uranus() {
        let orbit_normal = URANUS.orbit.normal();
        let north = URANUS.rotation_axis.to_direction();
        println!("orbit_normal: {}", orbit_normal);
        println!("north: {}", north);
        let expected = URANUS.axis_tilt;
        let actual = orbit_normal.angle_to(&north);
        println!("expected: {}, actual: {}", expected, actual);
        assert!(angle_eq_within(actual, expected, TILT_ACCURACY));
    }

    #[test]
    fn axis_tilt_of_neptune() {
        let orbit_normal = NEPTUNE.orbit.normal();
        let north = NEPTUNE.rotation_axis.to_direction();
        println!("orbit_normal: {}", orbit_normal);
        println!("north: {}", north);
        let expected = NEPTUNE.axis_tilt;
        let actual = orbit_normal.angle_to(&north);
        println!("expected: {}, actual: {}", expected, actual);
        assert!(angle_eq_within(actual, expected, TILT_ACCURACY));
    }

    #[test]
    fn axis_tilt_of_pluto() {
        let orbit_normal = PLUTO.orbit.normal();
        let north = PLUTO.rotation_axis.to_direction();
        println!("orbit_normal: {}", orbit_normal);
        println!("north: {}", north);
        let expected = PLUTO.axis_tilt;
        let actual = orbit_normal.angle_to(&north);
        println!("expected: {}, actual: {}", expected, actual);
        assert!(angle_eq_within(actual, expected, TILT_ACCURACY));
    }
}
