use crate::*;
use rand::Rng;

pub struct OrbitGenerator {
    pub stellar_mass: Mass,
    pub semi_major_axis: Length,
}

impl OrbitGenerator {
    pub fn gen_next<R: Rng>(&mut self, rng: &mut R) -> EllipticalOrbit {
        let semi_major_axis = self.semi_major_axis;
        self.semi_major_axis *= rng.gen_range(1.3..2.0);

        EllipticalOrbit::new(
            self.stellar_mass,
            semi_major_axis,
            rng.gen(),
            rng.gen(),
            rng.gen(),
        )
    }
}
