use physics_types::{Angle, Distance, Duration, Length, Mass, Squared, TimeFloat};
use std::ops::{Add, Mul, Sub};

// TODO orbit positions to be refactored. Calculate orbit path rather than doing the full Keplerian calculation each step.
// graph in Excel and consider the best way to approximate the calculation from precalculated values.

pub const GRAVITY_CONST: f64 = 6.674_08e-11;
pub const GRAVITY_CONST_SQRT: f64 = 8.169_504_268_926e-6;

#[derive(Debug, Default, Copy, Clone, PartialOrd, PartialEq)]
pub struct Eccentricity(f64);

impl Eccentricity {
    #[inline]
    pub fn new(value: f64) -> Self {
        // I don't know whether my equations work for high eccentricities
        if value < 0.0 || value > 0.9 {
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

#[derive(Debug, Copy, Clone)]
pub struct KeplerOrbit {
    pub period: Duration,
    pub semi_major_axis: Length,
    pub eccentricity: Eccentricity,
    pub eccentricity_angle: Angle,
    pub anomaly_offset: MeanAnomaly,
}

impl KeplerOrbit {
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
    pub fn position(&self, time: TimeFloat) -> Distance {
        if self.eccentricity.0 == 0.0 {
            circular_orbit_position(
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
    pub fn position(&self, time: TimeFloat) -> Distance {
        circular_orbit_position(time, self.period, self.radius, self.offset)
    }
}

#[inline]
fn circular_orbit_position(
    time: TimeFloat,
    period: Duration,
    radius: Length,
    offset: Angle,
) -> Distance {
    let angle = Angle::TAU * (time / period) + offset;
    Distance::from_angle_and_radius(angle, radius)
}

//
// impl World {
//
//     pub fn draw_body_orbit(
//         &self,
//         mesh_builder: &mut MeshBuilder,
//         body: Id<Body>,
//         color: Color,
//         camera: &crate::cameras::Camera,
//     ) -> Option<()> {
//         let orbit = self.state.body.orbit.get_opt(body)?;
//         self.draw_orbit(mesh_builder, *orbit, color, camera)
//     }
//
//     pub fn draw_orbit(
//         &self,
//         mesh_builder: &mut MeshBuilder,
//         orbit: Id<Orbit>,
//         color: Color,
//         camera: &crate::cameras::Camera,
//     ) -> Option<()> {
//         let system_view = &self.state.system_view;
//         let is_on_screen = |p: &&Point| system_view.is_on_screen_point(**p);
//
//         let points = self.calculate_orbit_points(orbit, camera);
//
//         points.iter()
//             .zip(points.iter().skip(1))
//             .enumerate()
//             .filter(|(_, (a, b))| is_on_screen(a) || is_on_screen(b))
//             .for_each(|(i, (a, b))| {
//                 let brightness = 1.0 - (i as f32 / TRAIL_POINT_COUNT as f32);
//                 let points = &[*a, *b];
//                 let color = Color::new(
//                     color.r * brightness,
//                     color.g * brightness,
//                     color.b * brightness,
//                     color.a,
//                 );
//                 let _ = mesh_builder.line(points, 0.5, color);
//             });
//
//         Some(())
//     }
//
//     fn calculate_orbit_points(
//         &self,
//         orbit: Id<Orbit>,
//         camera: &crate::cameras::Camera,
//     ) -> Vec<Point> {
//         let system_view = &self.state.system_view;
//         let orbit_params = &self.state.orbit.parameters[orbit];
//         let parent_orbit = self.state.orbit.parent[orbit];
//
//         let parent_pos = parent_orbit
//             .map(|p| self.calculate_absolute_position(p))
//             .unwrap_or_else(Position::zero);
//
//         let end_time = self.state.time.get_game_time() + orbit_params.period * RENDER_ORBIT_FRACTION.ceil();
//         let time_step = orbit_params.period / TRAIL_POINT_COUNT as f64;
//
//         (0..TRAIL_POINT_COUNT)
//             .map(|i| {
//                 let t = end_time - time_step * i as f64 * RENDER_ORBIT_FRACTION;
//                 let p = calculate_relative_position(orbit_params, t) + parent_pos;
//                 let px = camera.get_pixel_from_position(p, system_view);
//                 get_point(px)
//             })
//             .collect()
//     }
//
//     pub fn draw_if_child_orbit(
//         &self,
//         mesh_builder: &mut MeshBuilder,
//         body: Id<Body>,
//         color: Color,
//         camera: &Camera,
//         parent_orbit: Id<Orbit>,
//     ) -> Option<()> {
//         let orbit = self.state.body.orbit.get_opt(body)?;
//         let orbit_parent = self.state.orbit.parent.get_opt(orbit)?;
//
//         if *orbit_parent == parent_orbit {
//             self.draw_body_orbit(mesh_builder, body, color, camera);
//         }
//
//         Some(())
//     }
//
//     fn calculate_absolute_position(&self, orbit: Id<Orbit>) -> Position {
//         calculate_absolute_position(
//             orbit,
//             &self.state.orbit.relative_pos,
//             &self.state.orbit.parent,
//         )
//     }
//
//     pub fn calculate_absolute_velocity(&self, orbit: Id<Orbit>) -> Velocity {
//         self.calculate_absolute_velocity_at_time(orbit, self.state.time.get_game_time())
//     }
//
//     pub fn calculate_absolute_velocity_at_time(&self, orbit: Id<Orbit>, time: Time) -> Velocity {
//         calculate_absolute_velocity(
//             time,
//             orbit,
//             &self.state.orbit.parameters,
//             &self.state.orbit.parent,
//         )
//     }
// }
//
// fn calculate_absolute_position(
//     orbit: Id<Orbit>,
//     orbit_positions: &Component<Orbit, Position>,
//     orbit_parent: &Component<Orbit, Option<Id<Orbit>>>,
// ) -> Position {
//     let mut position = Position::zero();
//     let mut orbit = Some(orbit);
//
//     while let Some(o) = orbit {
//         if let Some((rel_pos, parent)) = (orbit_positions, orbit_parent).get(o) {
//             position += *rel_pos;
//             orbit = *parent;
//         }
//     }
//
//     position
// }
//
// fn calculate_relative_velocity(params: &OrbitParameters, time: Time) -> Velocity {
//     let dt = Time::in_seconds(1.0);
//     let p1 = calculate_relative_position(params, time);
//     let p2 = calculate_relative_position(params, time + dt);
//     (p2 - p1) / dt
// }
//
// fn calculate_absolute_velocity(
//     time: Time,
//     orbit: Id<Orbit>,
//     orbit_params: &Component<Orbit, OrbitParameters>,
//     orbit_parent: &Component<Orbit, Option<Id<Orbit>>>,
// ) -> Velocity {
//     let mut velocity = Velocity::zero();
//     let mut orbit = Some(orbit);
//
//     while let Some(o) = orbit {
//         if let Some((params, parent)) = (orbit_params, orbit_parent).get(o) {
//             velocity += calculate_relative_velocity(params, time);
//             orbit = *parent;
//         }
//     }
//
//     velocity
// }

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
        assert_eq!(40, std::mem::size_of::<KeplerOrbit>());
        assert_eq!(24, std::mem::size_of::<CircularOrbit>());
    }
}
