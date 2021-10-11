use crate::calc::*;
use physics_types::{Angle, Distance, Duration, Length, Mass, Polar, Squared, TimeFloat};

pub mod calc;
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

    pub fn circular_from_parent(parent: Mass, radius: Length, offset: Angle) -> Self {
        let period = orbital_period(radius, parent);
        Self::circular(period, radius, offset)
    }

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

impl std::ops::Mul<f64> for Eccentricity {
    type Output = f64;

    #[inline]
    fn mul(self, rhs: f64) -> Self::Output {
        self.0 * rhs
    }
}
