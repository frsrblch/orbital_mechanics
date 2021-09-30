use physics_types::{Angle, Length};
use std::ops::{Add, Sub};

/// The angle φ is in the range [0..π], and represents the angle relative to the poles
#[derive(Debug, Default, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct Phi(Angle);

impl From<f64> for Phi {
    fn from(fraction: f64) -> Self {
        Self(Angle::acos(1.0 - 2.0 * fraction))
    }
}

/// The angle θ represents the rotation of the spiral in the interval [0..Rτ]
/// Where R is the number of rotations, as calculated from the number of nodes by the `rotations` function
#[derive(Debug, Default, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct Theta(Angle);

impl Theta {
    pub(crate) fn fraction(fraction: f64, rotations: f64) -> Self {
        Self::rotations(Phi::from(fraction), rotations)
    }

    pub(crate) fn rotations(phi: Phi, rotations: f64) -> Self {
        Self(phi.0 * rotations)
    }
}

impl Add<Angle> for Theta {
    type Output = Theta;

    fn add(self, rhs: Angle) -> Self::Output {
        Theta(self.0 + rhs)
    }
}

impl Sub<Angle> for Theta {
    type Output = Theta;

    fn sub(self, rhs: Angle) -> Self::Output {
        Theta(self.0 - rhs)
    }
}

/// Represents a point on a sphere of arbitrary radius
#[derive(Debug, Default, Copy, Clone, PartialOrd, PartialEq)]
pub struct SphericalCoordinate {
    pub phi: Phi,
    pub theta: Theta,
}

impl SphericalCoordinate {
    pub fn position(&self, rho: Length) -> (Length, Length, Length) {
        let x = rho * self.phi.0.sin() * self.theta.0.cos();
        let y = rho * self.phi.0.sin() * self.theta.0.sin();
        let z = rho * self.phi.0.cos();
        (x, y, z)
    }

    pub fn unit_vector(&self) -> (f64, f64, f64) {
        let x = self.phi.0.sin() * self.theta.0.cos();
        let y = self.phi.0.sin() * self.theta.0.sin();
        let z = self.phi.0.cos();
        (x, y, z)
    }
}
