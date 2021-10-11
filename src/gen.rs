use crate::*;
use rand::distributions::uniform::SampleRange;
use rand::Rng;
use std::f64::consts::PI;

pub struct OrbitGenerator<R> {
    pub stellar_mass: Mass,
    pub semi_major_axis: Length,
    pub rng: R,
}

impl<R: Rng> Iterator for OrbitGenerator<R> {
    type Item = EllipticalOrbit;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.create_orbit())
    }
}

impl<R: Rng> OrbitGenerator<R> {
    pub fn create_orbit(&mut self) -> EllipticalOrbit {
        let semi_major_axis = self.semi_major_axis;
        self.advance_to_next();

        let eccentricity = self.eccentricity();
        let angle = Angle::in_rad(self.rng.gen_range(-PI..PI));
        let offset = Angle::in_rad(self.rng.gen_range(-PI..PI));

        EllipticalOrbit::new(
            self.stellar_mass,
            semi_major_axis,
            eccentricity,
            angle,
            offset,
        )
    }

    pub fn advance_to_next(&mut self) {
        self.semi_major_axis *= self.rng.gen_range(1.3..2.0);
    }

    fn eccentricity(&mut self) -> Eccentricity {
        let min = -0.1_f64.powf(1. / 3.);
        let max = 0.25_f64.powf(1. / 3.);
        let range = min..max;
        let e = range.sample_single(&mut self.rng).abs().powi(3);
        Eccentricity::new(e)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::thread_rng;

    #[test]
    fn gen_eccentricity() {
        let mut p = OrbitGenerator {
            stellar_mass: Mass::in_kg(1.0),
            semi_major_axis: Length::in_m(1.0),
            rng: thread_rng(),
        };

        for _ in 0..100 {
            let _e = p.eccentricity();
        }
    }
}
