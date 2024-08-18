pub trait PassiveRotation<T> {
    fn passive_rotation_to_new_z_axis(&self, new_z: &T) -> T;
}
