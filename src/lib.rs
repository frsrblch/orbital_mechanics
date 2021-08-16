use physics_types::{Angle, Distance, Duration, Length, Mass, Squared, TimeFloat};
use std::ops::{Add, Mul, Sub};

// TODO orbit positions to be refactored. Calculate orbit path rather than doing the full Keplerian calculation each step.
// graph in Excel and consider the best way to approximate the calculation from precalculated values.

pub const GRAVITY_CONST: f64 = 6.674_08e-11;
pub const GRAVITY_CONST_SQRT: f64 = 8.169_504_268_926e-6;

#[derive(Debug, Copy, Clone)]
pub struct EllipticalOrbit {
    pub period: Duration,
    pub semi_major_axis: Length,
    pub eccentricity: Eccentricity,
    pub eccentricity_angle: Angle,
    pub anomaly_offset: MeanAnomaly,
}

impl EllipticalOrbit {
    pub fn circular(period: Duration, radius: Length, offset: Angle) -> Self {
        Self {
            period,
            semi_major_axis: radius,
            eccentricity: Default::default(),
            eccentricity_angle: Default::default(),
            anomaly_offset: MeanAnomaly(offset),
        }
    }

    #[inline]
    pub fn distance(&self, time: TimeFloat) -> Distance {
        if self.eccentricity.0 == 0.0 {
            circular_orbit_distance(
                time,
                self.period,
                self.semi_major_axis,
                self.anomaly_offset.0,
            )
        } else {
            let mean_anomaly = MeanAnomaly::calculate(self.anomaly_offset, self.period, time);
            let eccentric_anomaly = EccentricAnomaly::calculate(mean_anomaly, self.eccentricity);
            let true_anomaly = TrueAnomaly::calculate(eccentric_anomaly, self.eccentricity);

            let radius = radius(true_anomaly, self.eccentricity, self.semi_major_axis);
            let angle = true_anomaly + self.eccentricity_angle;

            Distance::from_angle_and_radius(angle, radius)
        }
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
    pub fn distance(&self, time: TimeFloat) -> Distance {
        circular_orbit_distance(time, self.period, self.radius, self.offset)
    }
}

#[derive(Debug, Default, Copy, Clone, PartialOrd, PartialEq)]
pub struct Eccentricity(f64);

impl Eccentricity {
    #[inline]
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

impl Sub<f64> for Eccentricity {
    type Output = f64;

    #[inline]
    fn sub(self, rhs: f64) -> Self::Output {
        self.0 - rhs
    }
}

impl Sub<Eccentricity> for f64 {
    type Output = f64;

    #[inline]
    fn sub(self, rhs: Eccentricity) -> Self::Output {
        self - rhs.0
    }
}

impl Mul<f64> for Eccentricity {
    type Output = f64;

    #[inline]
    fn mul(self, rhs: f64) -> Self::Output {
        self.0 * rhs
    }
}

impl Mul<Eccentricity> for f64 {
    type Output = f64;

    #[inline]
    fn mul(self, rhs: Eccentricity) -> Self::Output {
        self * rhs.0
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq, PartialOrd)]
pub struct MeanAnomaly(pub Angle);

impl MeanAnomaly {
    #[inline]
    pub fn calculate(offset: MeanAnomaly, orbital_period: Duration, time: TimeFloat) -> Self {
        let orbit_fraction = time / orbital_period;
        let angle = orbit_fraction * Angle::TAU + offset;
        MeanAnomaly(angle)
    }

    #[inline]
    pub fn cos(self) -> f64 {
        self.0.cos()
    }

    #[inline]
    pub fn sin(self) -> f64 {
        self.0.sin()
    }
}

impl Add<Angle> for MeanAnomaly {
    type Output = Angle;

    #[inline]
    fn add(self, rhs: Angle) -> Self::Output {
        self.0 + rhs
    }
}

impl Add<MeanAnomaly> for Angle {
    type Output = Angle;

    #[inline]
    fn add(self, rhs: MeanAnomaly) -> Self::Output {
        self + rhs.0
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq, PartialOrd)]
pub struct EccentricAnomaly(pub Angle);

impl EccentricAnomaly {
    const ITERATION_COUNT: u8 = 5;

    #[inline]
    pub fn calculate(mean_anomaly: MeanAnomaly, eccentricity: Eccentricity) -> Self {
        // set E_0 to an appropriate initial value
        let mut ecc_anomaly = mean_anomaly
            + Angle::in_rad(
                eccentricity * mean_anomaly.sin()
                    / (1.0 - (mean_anomaly + Angle::in_rad(eccentricity.0)).sin()
                        + mean_anomaly.sin()),
            );

        // Newton-Raphson method
        for _ in 0..Self::ITERATION_COUNT {
            let f = ecc_anomaly.value - eccentricity * ecc_anomaly.sin() - mean_anomaly.0.value;
            let f_prime = 1.0 - eccentricity * ecc_anomaly.cos();
            ecc_anomaly -= Angle::in_rad(f / f_prime);
        }

        EccentricAnomaly(ecc_anomaly)
    }

    #[inline]
    pub fn cos(self) -> f64 {
        self.0.cos()
    }

    #[inline]
    pub fn sin(self) -> f64 {
        self.0.sin()
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq, PartialOrd)]
pub struct TrueAnomaly(pub Angle);

impl Add<Angle> for TrueAnomaly {
    type Output = Angle;

    #[inline]
    fn add(self, rhs: Angle) -> Self::Output {
        self.0 + rhs
    }
}

impl Add<TrueAnomaly> for Angle {
    type Output = Angle;

    #[inline]
    fn add(self, rhs: TrueAnomaly) -> Self::Output {
        self + rhs.0
    }
}

impl TrueAnomaly {
    #[inline]
    pub fn calculate(eccentric_anomaly: EccentricAnomaly, eccentricity: Eccentricity) -> Self {
        let cos = eccentric_anomaly.cos();
        let anomaly = Angle::in_rad(((cos - eccentricity) / (1.0 - eccentricity * cos)).acos());
        if eccentric_anomaly.0 % Angle::TAU > Angle::PI {
            TrueAnomaly(Angle::TAU - anomaly)
        } else {
            TrueAnomaly(anomaly)
        }
    }

    #[inline]
    pub fn cos(self) -> f64 {
        self.0.cos()
    }

    #[inline]
    pub fn sin(self) -> f64 {
        self.0.sin()
    }
}

#[inline]
pub fn radius(
    true_anomaly: TrueAnomaly,
    eccentricity: Eccentricity,
    semi_major_axis: Length,
) -> Length {
    semi_major_axis * (1.0 - eccentricity.squared()) / (1.0 + eccentricity * true_anomaly.cos())
}

pub fn orbital_period(semi_major_axis: Length, parent_mass: Mass) -> Duration {
    const MULTIPLIER: f64 = Angle::TAU.value / GRAVITY_CONST_SQRT;
    let value = MULTIPLIER * (semi_major_axis.value.powi(3) / parent_mass.value).sqrt();
    Duration::in_s(value)
}

#[inline]
fn circular_orbit_distance(
    time: TimeFloat,
    period: Duration,
    radius: Length,
    offset: Angle,
) -> Distance {
    let angle = Angle::TAU * (time / period) + offset;
    Distance::from_angle_and_radius(angle, radius)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn earth_orbital_period() {
        let earth_sma = Length::in_m(149598023e3);
        let m_sun = Mass::in_kg(1.9885e30);
        let period = orbital_period(earth_sma, m_sun);

        assert!(period > Duration::in_days(365.2));
        assert!(period < Duration::in_days(365.3));
    }

    #[test]
    fn orbit_sizes() {
        assert_eq!(40, std::mem::size_of::<EllipticalOrbit>());
        assert_eq!(24, std::mem::size_of::<CircularOrbit>());
    }
}
