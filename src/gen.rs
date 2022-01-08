use crate::*;
use rand::Rng;
use std::cell::Cell;

pub struct OrbitGenerator {
    pub stellar_mass: Mass,
    pub semi_major_axis: Cell<Length>,
}

impl Distribution<EllipticalOrbit> for OrbitGenerator {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> EllipticalOrbit {
        let semi_major_axis = self.semi_major_axis.get();
        self.semi_major_axis
            .set(self.semi_major_axis.get() * rng.gen_range(1.3..2.0));

        EllipticalOrbit::new(
            self.stellar_mass,
            semi_major_axis,
            rng.gen(),
            rng.gen(),
            rng.gen(),
        )
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::thread_rng;

    #[test]
    fn each_orbit_larger_than_the_last() {
        let generator = OrbitGenerator {
            stellar_mass: Mass::in_kg(1.0),
            semi_major_axis: Cell::new(Length::in_m(1.0)),
        };

        let rng = &mut thread_rng();

        let mut orbit = generator.sample(rng);

        for _ in 0..10 {
            let next = generator.sample(rng);
            assert!(next.semi_major_axis > orbit.semi_major_axis);
            orbit = next;
        }
    }
}
