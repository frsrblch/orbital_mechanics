use orbital_mechanics::{Eccentricity, EllipticalOrbit};
use physics_types::{Duration, Length, TimeFloat};
use std::io::Write;

fn main() {
    let orbit = EllipticalOrbit {
        period: Duration::in_s(1.0),
        semi_major_axis: Length::in_m(1.0),
        eccentricity: Eccentricity::new(0.4),
        eccentricity_angle: Default::default(),
        offset: Default::default(),
    };

    let mut file = std::fs::File::create("orbit.txt").unwrap();

    const N: usize = 1000;
    for n in 0..=N {
        let t = n as f64 / N as f64;
        let t = TimeFloat::in_s(t);
        let r = orbit.radius(t);
        writeln!(&mut file, "{}\t{}", t.value.value, r.value).unwrap();
    }
}
