//! Cartesian coordinates in 3D space.

use serde::{Deserialize, Serialize};
use simple_si_units::{
    base::Distance,
    geometry::{Angle, Area},
};
use std::{
    f64::consts::PI,
    fmt::Display,
    ops::{Add, Div, Mul, Neg, Sub},
};

use crate::{
    angle_helper::{safe_acos, HALF_CIRC},
    earth_equatorial::EarthEquatorial,
    equatorial::Equatorial,
    error::AstroCoordsError,
    NORMALIZATION_THRESHOLD,
};

use super::{
    direction::Direction, ecliptic::Ecliptic, spherical::Spherical, transformations::rotations::*,
};

/// Cartesian coordinates in 3D space.
///
/// The members x, y and z are typed as `Distance<f64>`.
///
/// The struct implments basic arithmetic operations, such as addition, subtraction, multiplication and division.
///
/// # Examples
/// ```
/// use simple_si_units::base::Distance;
/// use astro_coords::cartesian::Cartesian;
///
/// let c1 = Cartesian::new(Distance::from_meters(1.), Distance::from_meters(2.), Distance::from_meters(3.));
/// let c2 = Cartesian::new(Distance::from_meters(4.), Distance::from_meters(5.), Distance::from_meters(6.));
///
/// let add = &c1 + &c2;
/// let expected = Cartesian::new(Distance::from_meters(5.), Distance::from_meters(7.), Distance::from_meters(9.));
/// assert!(add.eq_within(&expected, Distance::from_meters(1e-5)));
///
/// let sub = &c1 - &c2;
/// let expected = Cartesian::new(Distance::from_meters(-3.), Distance::from_meters(-3.), Distance::from_meters(-3.));
/// assert!(sub.eq_within(&expected, Distance::from_meters(1e-5)));
///
/// let mul = &c1 * 2.;
/// let expected = Cartesian::new(Distance::from_meters(2.), Distance::from_meters(4.), Distance::from_meters(6.));
/// assert!(mul.eq_within(&expected, Distance::from_meters(1e-5)));
///
/// let div = &c1 / 2.;
/// let expected = Cartesian::new(Distance::from_meters(0.5), Distance::from_meters(1.), Distance::from_meters(1.5));
/// assert!(div.eq_within(&expected, Distance::from_meters(1e-5)));
///
/// let neg = -c1;
/// let expected = Cartesian::new(Distance::from_meters(-1.), Distance::from_meters(-2.), Distance::from_meters(-3.));
/// assert!(neg.eq_within(&expected, Distance::from_meters(1e-5)));
/// ```
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Cartesian {
    /// The x-ordinate.
    pub x: Distance<f64>,
    /// The y-ordinate.
    pub y: Distance<f64>,
    /// The z-ordinate.
    pub z: Distance<f64>,
}

impl Cartesian {
    /// The origin of the coordinate system, given by (x,y,z)=(0,0,0).
    pub const ORIGIN: Cartesian = Cartesian {
        x: Distance { m: 0. },
        y: Distance { m: 0. },
        z: Distance { m: 0. },
    };

    /// Creates a new Cartesian object.
    pub const fn new(x: Distance<f64>, y: Distance<f64>, z: Distance<f64>) -> Cartesian {
        Cartesian { x, y, z }
    }

    /// Tests if the coordinates are equal to another set of coordinates within a certain accuracy.
    ///
    /// # Examples
    /// ```
    /// use simple_si_units::base::Distance;
    /// use astro_coords::cartesian::Cartesian;
    ///
    /// let c1 = Cartesian::new(Distance::from_meters(1.), Distance::from_meters(2.), Distance::from_meters(3.));
    /// let c2 = Cartesian::new(Distance::from_meters(1.), Distance::from_meters(2.), Distance::from_meters(3.));
    /// assert!(c1.eq_within(&c2, Distance::from_meters(0.1)));
    /// ```
    pub fn eq_within(&self, other: &Cartesian, accuracy: Distance<f64>) -> bool {
        (self.x.m - other.x.m).abs() < accuracy.m
            && (self.y.m - other.y.m).abs() < accuracy.m
            && (self.z.m - other.z.m).abs() < accuracy.m
    }

    /// Returns the length of the vector defined by the coordinates.
    ///
    /// Note that this operation includes a square root calculation, which is computationally expensive.
    /// If performance is an issue and you only need to compare the lengths of two vectors, consider using `length_squared()` instead.
    ///
    /// # Examples
    /// ```
    /// use simple_si_units::base::Distance;
    /// use astro_coords::cartesian::Cartesian;
    ///
    /// let coordinates = Cartesian::new(Distance::from_meters(3.), Distance::from_meters(4.), Distance::from_meters(5.));
    /// assert!((coordinates.length().to_m() - 7.0710678118654755).abs() < 1e-5);
    /// ```
    pub fn length(&self) -> Distance<f64> {
        let x = self.x.m;
        let y = self.y.m;
        let z = self.z.m;
        Distance {
            m: (x * x + y * y + z * z).sqrt(),
        }
    }

    /// Returns the square of the length of the vector defined by the coordinates.
    ///
    /// # Examples
    /// ```
    /// use simple_si_units::base::Distance;
    /// use astro_coords::cartesian::Cartesian;
    ///
    /// let coordinates = Cartesian::new(Distance::from_meters(3.), Distance::from_meters(4.), Distance::from_meters(5.));
    /// assert!((coordinates.length_squared().to_m2() - 50.).abs() < 1e-5);
    /// ```
    pub fn length_squared(&self) -> Area<f64> {
        let x = self.x.m;
        let y = self.y.m;
        let z = self.z.m;
        Area {
            m2: x * x + y * y + z * z,
        }
    }

    /// Returns the distance between two sets of Cartesian coordinates.
    ///
    /// # Examples
    /// ```
    /// use simple_si_units::base::Distance;
    /// use astro_coords::cartesian::Cartesian;
    ///
    /// let c1 = Cartesian::new(Distance::from_meters(1.), Distance::from_meters(2.), Distance::from_meters(3.));
    /// let c2 = Cartesian::new(Distance::from_meters(4.), Distance::from_meters(5.), Distance::from_meters(6.));
    /// assert!((c1.distance(&c2).to_m() - 5.196152422706632).abs() < 1e-5);
    /// ```
    pub fn distance(&self, other: &Cartesian) -> Distance<f64> {
        let diff = self - other;
        diff.length()
    }

    /// Returns cartesian coordinates rotated around an axis by a certain angle.
    ///
    /// # Examples
    /// ```
    /// use simple_si_units::base::Distance;
    /// use simple_si_units::geometry::Angle;
    /// use astro_coords::cartesian::Cartesian;
    /// use astro_coords::direction::Direction;
    ///
    /// let coordinates = Cartesian::new(Distance::from_meters(1.), Distance::from_meters(0.), Distance::from_meters(0.));
    /// let angle = Angle::from_degrees(90.);
    /// let axis = Direction::new(0., 0., 1.).unwrap();
    /// let rotated = coordinates.rotated(angle, &axis);
    /// let expected = Cartesian::new(Distance::from_meters(0.), Distance::from_meters(1.), Distance::from_meters(0.));
    /// assert!(rotated.eq_within(&expected, Distance::from_meters(1e-5)));
    /// ```
    pub fn rotated(&self, angle: Angle<f64>, axis: &Direction) -> Cartesian {
        let (x, y, z) = rotated_tuple((self.x, self.y, self.z), angle, axis);
        Cartesian { x, y, z }
    }

    pub fn rotated_x(&self, angle: Angle<f64>) -> Cartesian {
        let (x, y, z) = rotated_x_tuple((self.x, self.y, self.z), angle);
        Cartesian { x, y, z }
    }

    pub fn rotated_y(&self, angle: Angle<f64>) -> Cartesian {
        let (x, y, z) = rotated_y_tuple((self.x, self.y, self.z), angle);
        Cartesian { x, y, z }
    }

    pub fn rotated_z(&self, angle: Angle<f64>) -> Cartesian {
        let (x, y, z) = rotated_z_tuple((self.x, self.y, self.z), angle);
        Cartesian { x, y, z }
    }

    /// Returns the angle between two sets of Cartesian coordinates.
    ///
    /// # Examples
    /// ```
    /// use simple_si_units::base::Distance;
    /// use simple_si_units::geometry::Angle;
    /// use astro_coords::cartesian::Cartesian;
    ///
    /// let c1 = Cartesian::new(Distance::from_meters(1.), Distance::from_meters(0.), Distance::from_meters(0.));
    /// let c2 = Cartesian::new(Distance::from_meters(0.), Distance::from_meters(1.), Distance::from_meters(0.));
    /// let angle = c1.angle_to(&c2).unwrap();
    /// assert!((angle.to_degrees() - 90.).abs() < 1e-5);
    /// ```
    pub fn angle_to(&self, other: &Cartesian) -> Result<Angle<f64>, AstroCoordsError> {
        // Use that
        // cos²(θ) = (a · b)² / (|a|² |b|²) = 1/2 + cos(2 θ) / 2
        // => cos(2 θ) = 2 (a · b)² / (|a|² |b|²) - 1
        // => θ = 1/2 arccos(2 (a · b)² / (|a|² |b|²) - 1)
        let len_sqrd_a = self.length_squared().m2;
        let len_sqrd_b = other.length_squared().m2;
        let len_sqrd_prod = len_sqrd_a * len_sqrd_b;
        if len_sqrd_prod < NORMALIZATION_THRESHOLD {
            return Err(AstroCoordsError::NormalizingZeroVector);
        }
        let dot_prod = self.x.m * other.x.m + self.y.m * other.y.m + self.z.m * other.z.m;
        let cos_2_theta = 2. * dot_prod * dot_prod / len_sqrd_prod - 1.;
        let angle = safe_acos(cos_2_theta) / 2.;
        if dot_prod >= 0. {
            Ok(angle)
        } else {
            Ok(HALF_CIRC - angle)
        }
    }

    /// Returns a unit vector in the direction of the coordinates.
    ///
    /// For coordinates with a very small length, this function will return an error.
    ///
    /// # Examples
    /// ```
    /// use simple_si_units::base::Distance;
    /// use astro_coords::cartesian::Cartesian;
    ///
    /// let coordinates = Cartesian::new(Distance::from_meters(1.), Distance::from_meters(1.), Distance::from_meters(0.));
    /// let direction = coordinates.to_direction().unwrap();
    /// assert!((direction.x() - 0.7071067811865476).abs() < 1e-5);
    /// assert!((direction.y() - 0.7071067811865476).abs() < 1e-5);
    /// assert!((direction.z() - 0.).abs() < 1e-5);
    /// ```
    pub fn to_direction(&self) -> Result<Direction, AstroCoordsError> {
        Direction::new(self.x.m, self.y.m, self.z.m)
    }

    pub fn to_earth_equatorial(&self) -> Result<EarthEquatorial, AstroCoordsError> {
        Ok(self.to_direction()?.to_earth_equatorial())
    }

    pub fn to_ecliptic(&self) -> Result<Ecliptic, AstroCoordsError> {
        Ok(Ecliptic {
            spherical: self.to_spherical()?,
        })
    }

    pub fn to_equatorial(&self, rotation_axis: Direction) -> Result<Equatorial, AstroCoordsError> {
        let spherical = self
            .to_direction()?
            .passive_rotation_to_new_z_axis(&rotation_axis)
            .to_spherical();
        Ok(Equatorial::new(spherical, rotation_axis))
    }

    /// Returns spherical coordinates representing the direction of the coordinates.
    ///
    /// For coordinates with a very small length, this function will return an error.
    ///
    /// # Examples
    /// ```
    /// use simple_si_units::base::Distance;
    /// use astro_coords::cartesian::Cartesian;
    ///
    /// let coordinates = Cartesian::new(Distance::from_meters(1.), Distance::from_meters(1.), Distance::from_meters(0.));
    /// let spherical = coordinates.to_spherical().unwrap();
    /// assert!((spherical.longitude.to_degrees() - 45.).abs() < 1e-5);
    /// assert!((spherical.latitude.to_degrees() - 0.).abs() < 1e-5);
    /// ```
    pub fn to_spherical(&self) -> Result<Spherical, AstroCoordsError> {
        Spherical::cartesian_to_spherical((self.x.m, self.y.m, self.z.m))
    }
}

impl Add for &Cartesian {
    type Output = Cartesian;

    fn add(self, other: &Cartesian) -> Cartesian {
        Cartesian {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Add for Cartesian {
    type Output = Cartesian;

    fn add(self, other: Cartesian) -> Cartesian {
        Cartesian {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for &Cartesian {
    type Output = Cartesian;

    fn sub(self, other: &Cartesian) -> Cartesian {
        Cartesian {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Sub for Cartesian {
    type Output = Cartesian;

    fn sub(self, other: Cartesian) -> Cartesian {
        Cartesian {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Mul<f64> for &Cartesian {
    type Output = Cartesian;

    fn mul(self, f: f64) -> Cartesian {
        Cartesian {
            x: self.x * f,
            y: self.y * f,
            z: self.z * f,
        }
    }
}

impl Mul<f64> for Cartesian {
    type Output = Cartesian;

    fn mul(self, f: f64) -> Cartesian {
        Cartesian {
            x: self.x * f,
            y: self.y * f,
            z: self.z * f,
        }
    }
}

impl Div<f64> for &Cartesian {
    type Output = Cartesian;

    fn div(self, f: f64) -> Cartesian {
        Cartesian {
            x: self.x / f,
            y: self.y / f,
            z: self.z / f,
        }
    }
}

impl Div<f64> for Cartesian {
    type Output = Cartesian;

    fn div(self, f: f64) -> Cartesian {
        Cartesian {
            x: self.x / f,
            y: self.y / f,
            z: self.z / f,
        }
    }
}

impl Neg for Cartesian {
    type Output = Cartesian;

    fn neg(self) -> Cartesian {
        Cartesian {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Neg for &Cartesian {
    type Output = Cartesian;

    fn neg(self) -> Cartesian {
        Cartesian {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Display for Cartesian {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let length = self.length();
        if length.to_parsec() > 0.1 {
            write!(
                f,
                "({:.2} pc, {:.2} pc, {:.2} pc)",
                self.x.to_parsec(),
                self.y.to_parsec(),
                self.z.to_parsec()
            )
        } else if length.to_au() > 0.1 {
            write!(
                f,
                "({:.2} AU, {:.2} AU, {:.2} AU)",
                self.x.to_au(),
                self.y.to_au(),
                self.z.to_au()
            )
        } else {
            write!(
                f,
                "({:.2} km, {:.2} km, {:.2} km)",
                self.x.to_km(),
                self.y.to_km(),
                self.z.to_km()
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_ACCURACY: f64 = 1e-5;

    #[test]
    fn test_length() {
        let coordinates = Cartesian {
            x: Distance::from_meters(3.),
            y: Distance::from_meters(4.),
            z: Distance::from_meters(5.),
        };

        assert!((coordinates.length().to_m() - 7.0710678118654755).abs() < TEST_ACCURACY);
    }

    #[test]
    fn subtraction_tests() {
        let ordinates = vec![0., 1., -1., 12.17];
        for x1 in ordinates.iter() {
            for y1 in ordinates.iter() {
                for z1 in ordinates.iter() {
                    for x2 in ordinates.iter() {
                        for y2 in ordinates.iter() {
                            for z2 in ordinates.iter() {
                                let c1 = Cartesian::new(
                                    Distance::from_meters(*x1),
                                    Distance::from_meters(*y1),
                                    Distance::from_meters(*z1),
                                );
                                let c2 = Cartesian::new(
                                    Distance::from_meters(*x2),
                                    Distance::from_meters(*y2),
                                    Distance::from_meters(*z2),
                                );
                                let c3 = &c1 - &c2;
                                let c4 = c1.clone() - c2.clone();
                                let c5 = c1 + -c2;
                                assert!(c3.eq_within(&c4, Distance { m: TEST_ACCURACY }));
                                assert!(c4.eq_within(&c5, Distance { m: TEST_ACCURACY }));
                            }
                        }
                    }
                }
            }
        }
    }
}
