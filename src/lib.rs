use physics_types::{Angle, Distance, Duration, Length, Mass, Polar, Squared, TimeFloat};
use std::ops::Mul;

// TODO orbit positions to be refactored. Calculate orbit path rather than doing the full Keplerian calculation each step.
// graph in Excel and consider the best way to approximate the calculation from precalculated values.

pub mod coordinates;
pub mod gen;

#[derive(Debug, Copy, Clone)]
pub struct EllipticalOrbit {
    pub period: Duration,
    pub semi_major_axis: Length,
    pub eccentricity: Eccentricity,
    pub eccentricity_angle: Angle,
    pub offset: Angle,
}

impl EllipticalOrbit {
    #[inline]
    pub fn new(
        parent: Mass,
        semi_major_axis: Length,
        eccentricity: Eccentricity,
        eccentricity_angle: Angle,
        offset: Angle,
    ) -> Self {
        Self {
            period: orbital_period(semi_major_axis, parent),
            semi_major_axis,
            eccentricity,
            eccentricity_angle,
            offset,
        }
    }

    #[inline]
    pub fn circular_from_parent(parent: Mass, radius: Length, offset: Angle) -> Self {
        let period = orbital_period(radius, parent);
        Self::circular(period, radius, offset)
    }

    #[inline]
    pub fn circular(period: Duration, radius: Length, offset: Angle) -> Self {
        Self {
            period,
            semi_major_axis: radius,
            eccentricity: Default::default(),
            eccentricity_angle: Default::default(),
            offset,
        }
    }

    #[inline]
    fn is_circular(&self) -> bool {
        self.eccentricity == Eccentricity::default()
    }

    #[inline]
    pub fn distance(&self, time: TimeFloat) -> Distance {
        self.polar(time).euclidean()
    }

    #[inline]
    pub fn polar(&self, time: TimeFloat) -> Polar<Length> {
        if self.is_circular() {
            let angle = circular_orbit_angle(time, self.period, self.offset);
            let radius = self.semi_major_axis;

            Polar {
                x: radius,
                y: angle,
            }
        } else {
            let true_anomaly = self.true_anomaly(time);

            let angle = true_anomaly.0 + self.eccentricity_angle;
            let radius = self.elliptical_radius(true_anomaly);

            Polar {
                x: radius,
                y: angle,
            }
        }
    }

    #[inline]
    pub fn radius(&self, time: TimeFloat) -> Length {
        if self.is_circular() {
            self.semi_major_axis
        } else {
            let true_anomaly = self.true_anomaly(time);

            self.elliptical_radius(true_anomaly)
        }
    }

    #[inline]
    pub fn angle(&self, time: TimeFloat) -> Angle {
        if self.is_circular() {
            circular_orbit_angle(time, self.period, self.offset)
        } else {
            let true_anomaly = self.true_anomaly(time);
            true_anomaly.0 + self.eccentricity_angle
        }
    }

    #[inline]
    fn true_anomaly(&self, time: TimeFloat) -> TrueAnomaly {
        let mean_anomaly = MeanAnomaly::calculate(self.offset, self.period, time);
        let eccentric_anomaly = EccentricAnomaly::calculate(mean_anomaly, self.eccentricity);
        TrueAnomaly::calculate(eccentric_anomaly, self.eccentricity)
    }

    #[inline]
    fn elliptical_radius(&self, true_anomaly: TrueAnomaly) -> Length {
        self.semi_major_axis * (1.0 - self.eccentricity.squared())
            / (1.0 + self.eccentricity * true_anomaly.0.cos())
    }
}

#[derive(Debug, Copy, Clone)]
pub struct CircularOrbit {
    period: Duration,
    radius: Length,
    offset: Angle,
}

impl CircularOrbit {
    #[inline]
    pub fn from_parent(parent: Mass, radius: Length, offset: Angle) -> Self {
        Self {
            period: orbital_period(radius, parent),
            radius,
            offset,
        }
    }

    #[inline]
    pub fn distance(&self, time: TimeFloat) -> Distance {
        circular_orbit_distance(time, self.period, self.radius, self.offset)
    }
}

#[derive(Debug, Default, Copy, Clone, PartialOrd, PartialEq)]
pub struct Eccentricity(f64);

impl Eccentricity {
    pub fn new(value: f64) -> Self {
        // I don't know whether these equations work for high eccentricities
        if !(0.0..0.9).contains(&value) {
            panic!("Invalid eccentricity: {}", value);
        }

        Eccentricity(value)
    }
}

impl Squared for Eccentricity {
    type Output = f64;

    #[inline]
    fn squared(self) -> Self::Output {
        self.0 * self.0
    }
}

impl Mul<f64> for Eccentricity {
    type Output = f64;

    #[inline]
    fn mul(self, rhs: f64) -> Self::Output {
        self.0 * rhs
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq, PartialOrd)]
struct MeanAnomaly(Angle);

impl MeanAnomaly {
    fn calculate(offset: Angle, orbital_period: Duration, time: TimeFloat) -> Self {
        let orbit_fraction = time / orbital_period;
        let angle = orbit_fraction * Angle::TAU + offset;
        MeanAnomaly(angle)
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq, PartialOrd)]
struct EccentricAnomaly(Angle);

impl EccentricAnomaly {
    const ITERATION_COUNT: u8 = 3;

    fn calculate(mean_anomaly: MeanAnomaly, eccentricity: Eccentricity) -> Self {
        // set E_0 to an appropriate initial value
        let ma_sin = mean_anomaly.0.sin();

        let mut ecc_anomaly = mean_anomaly.0
            + Angle::in_rad(
                eccentricity * ma_sin
                    / (1.0 - (mean_anomaly.0 + Angle::in_rad(eccentricity.0)).sin() + ma_sin),
            );

        // Newton-Raphson method
        for _ in 0..Self::ITERATION_COUNT {
            let (sin, cos) = ecc_anomaly.value.sin_cos();
            let f = ecc_anomaly.value - eccentricity * sin - mean_anomaly.0.value;
            let f_prime = 1.0 - eccentricity * cos;
            ecc_anomaly -= Angle::in_rad(f / f_prime);
        }

        EccentricAnomaly(ecc_anomaly)
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq, PartialOrd)]
struct TrueAnomaly(Angle);

impl TrueAnomaly {
    fn calculate(eccentric_anomaly: EccentricAnomaly, eccentricity: Eccentricity) -> Self {
        let cos = eccentric_anomaly.0.cos();
        let anomaly = Angle::in_rad(((cos - eccentricity.0) / (1.0 - eccentricity * cos)).acos());
        if eccentric_anomaly.0 % Angle::TAU > Angle::PI {
            TrueAnomaly(Angle::TAU - anomaly)
        } else {
            TrueAnomaly(anomaly)
        }
    }
}

#[inline]
pub fn orbital_period(semi_major_axis: Length, parent_mass: Mass) -> Duration {
    const GRAVITY_CONST_SQRT: f64 = 8.169_504_268_926e-6;
    const MULTIPLIER: f64 = Angle::TAU.value / GRAVITY_CONST_SQRT;
    let value = MULTIPLIER * (semi_major_axis.value.powi(3) / parent_mass.value).sqrt();
    Duration::in_s(value)
}

fn circular_orbit_distance(
    time: TimeFloat,
    period: Duration,
    radius: Length,
    offset: Angle,
) -> Distance {
    let angle = circular_orbit_angle(time, period, offset);
    Distance::from_angle_and_magnitude(angle, radius)
}

#[inline]
fn circular_orbit_angle(time: TimeFloat, period: Duration, offset: Angle) -> Angle {
    Angle::TAU * (time / period) + offset
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn earth_orbital_period() {
        let earth_sma = Length::in_m(149598023e3);
        let m_sun = Mass::in_kg(1.9885e30);
        let period = orbital_period(earth_sma, m_sun);

        assert!(period > Duration::in_d(365.2));
        assert!(period < Duration::in_d(365.3));
    }

    #[test]
    fn orbit_sizes() {
        assert_eq!(40, std::mem::size_of::<EllipticalOrbit>());
        assert_eq!(24, std::mem::size_of::<CircularOrbit>());
    }
}
