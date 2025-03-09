use uom::si::f64::Angle;

use crate::direction::Direction;

pub trait ActiveRotation<T> {
    fn rotated(&self, angle: Angle, axis: &Direction) -> T;

    fn rotated_x(&self, angle: Angle) -> T;

    fn rotated_y(&self, angle: Angle) -> T;

    fn rotated_z(&self, angle: Angle) -> T;

    fn active_rotation_to_new_z_axis(&self, new_z: &T) -> T;
}
