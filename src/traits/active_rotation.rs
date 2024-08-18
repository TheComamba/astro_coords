use simple_si_units::geometry::Angle;

pub trait ActiveRotation<T> {
    fn rotated(&self, angle: Angle<f64>, axis: &T) -> T;

    fn rotated_x(&self, angle: Angle<f64>) -> T;

    fn rotated_y(&self, angle: Angle<f64>) -> T;

    fn rotated_z(&self, angle: Angle<f64>) -> T;

    fn active_rotation_to_new_z_axis(&self, new_z: &T) -> T;
}
