use criterion::*;
use orbital_mechanics::{Eccentricity, EllipticalOrbit};
use physics_types::*;

criterion_main! {
    calc_group,
}

criterion_group! {
    calc_group,
    calc,
}

fn calc(criterion: &mut Criterion) {
    let orbit = EllipticalOrbit {
        period: Duration::in_s(1.0),
        semi_major_axis: Length::in_m(1.0),
        eccentricity: Eccentricity::new(0.1),
        eccentricity_angle: Angle::zero(),
        offset: Angle::zero(),
    };

    let circle = EllipticalOrbit::circular(Duration::in_s(1.0), Length::in_m(1.0), Angle::zero());

    let t_a = TimeIndex::in_s(0.83);
    let t_b = TimeIndex::in_s(0.87);
    let t = TimeIndex::in_s(0.85);

    let dt_inv = 1.0 / (t_b - t_a);

    let p_a = orbit.polar(t_a);
    let p_b = orbit.polar(t_b);

    let f = (p_b - p_a) * dt_inv;

    let p_calc = orbit.polar(t);
    let d_calc = p_calc.euclidean();
    println!("calc: {}, {}", d_calc.x.value, d_calc.y.value);

    let p_approx = p_a + (p_b - p_a) * dt_inv * (t - t_a);
    let d_approx = p_approx.euclidean();
    println!("approx: {}, {}", d_approx.x.value, d_approx.y.value);

    println!("diff: {}", (d_calc - d_approx).magnitude().value);

    criterion
        .bench_function("calc", |bench| {
            bench.iter(|| {
                black_box(orbit.distance(t));
            })
        })
        .bench_function("circle", |bench| {
            bench.iter(|| black_box(circle.radius(t)));
        })
        .bench_function("approx", |bench| {
            bench.iter(|| {
                let p = p_a + (p_b - p_a) * dt_inv * (t - t_a);
                black_box(p.euclidean());
            })
        })
        .bench_function("approx simplified", |bench| {
            bench.iter(|| {
                let p = p_a + f * (t - t_a);
                black_box(p.euclidean());
            })
        });
}
