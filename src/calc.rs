use super::*;

#[derive(Debug, Default, Copy, Clone, PartialEq, PartialOrd)]
pub struct MeanAnomaly(pub Angle);

impl MeanAnomaly {
    pub fn calculate(offset: Angle, orbital_period: Duration, time: TimeIndex) -> Self {
        let orbit_fraction = time / orbital_period;
        let angle = orbit_fraction * Angle::TAU + offset;
        Self(angle)
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq, PartialOrd)]
pub struct EccentricAnomaly(pub Angle);

impl EccentricAnomaly {
    const ITERATION_COUNT: u8 = 3;

    pub fn calculate(mean_anomaly: MeanAnomaly, eccentricity: Eccentricity) -> Self {
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
            ecc_anomaly.value -= f / f_prime;
        }

        Self(ecc_anomaly)
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq, PartialOrd)]
pub struct TrueAnomaly(pub Angle);

impl TrueAnomaly {
    pub fn calculate(eccentric_anomaly: EccentricAnomaly, eccentricity: Eccentricity) -> Self {
        let cos = eccentric_anomaly.0.cos();
        let anomaly = Angle::in_rad(((cos - eccentricity.0) / (1.0 - eccentricity * cos)).acos());

        if eccentric_anomaly.0 % Angle::TAU > Angle::PI {
            Self(Angle::TAU - anomaly)
        } else {
            Self(anomaly)
        }
    }
}

#[inline]
pub fn orbital_period(semi_major_axis: Length, parent_mass: Mass) -> Duration {
    const MULTIPLIER: f64 = Angle::TAU.value / 8.169_504_268_926e-6; // TAU / sqrt(G)
    let value = MULTIPLIER * (semi_major_axis.value.powi(3) / parent_mass.value).sqrt();
    Duration::in_s(value)
}

#[inline]
pub fn circular_orbit_distance(
    time: TimeIndex,
    period: Duration,
    radius: Length,
    offset: Angle,
) -> Distance {
    let angle = circular_orbit_angle(time, period, offset);
    Distance::from_angle_and_magnitude(angle, radius)
}

#[inline]
pub fn circular_orbit_angle(time: TimeIndex, period: Duration, offset: Angle) -> Angle {
    Angle::TAU * (time / period) + offset
}

#[cfg(test)]
mod test {
    use super::*;
    use physics_types::{AU, KG, YR};

    #[test]
    fn earth_orbital_period() {
        let earth_sma = 1.0 * AU;
        let m_sun = 1.9885e30 * KG;

        let period = orbital_period(earth_sma, m_sun);

        assert!(period > YR * 0.999);
        assert!(period < YR * 1.001);
    }

    #[test]
    fn orbit_sizes() {
        assert_eq!(40, std::mem::size_of::<EllipticalOrbit>());
        assert_eq!(24, std::mem::size_of::<CircularOrbit>());
    }
}
