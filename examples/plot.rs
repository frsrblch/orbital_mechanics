use orbital_mechanics::calc::*;
use orbital_mechanics::Eccentricity;
use physics_types::*;
use plotters::prelude::*;

const OUT_FILE_NAME: &'static str = "3d-plot.svg";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let area = SVGBackend::new(OUT_FILE_NAME, (1024, 760)).into_drawing_area();

    area.fill(&WHITE)?;

    let e_axis = (0.0..0.8).step(0.1);
    let ma_axix = (0.0..1.0).step(0.25);

    let mut chart = ChartBuilder::on(&area)
        // .caption(format!("3D Plot Test"), ("sans", 20))
        .build_cartesian_3d(e_axis.clone(), 0.0..1.0, ma_axix.clone())?;

    chart.with_projection(|mut pb| {
        pb.yaw = 2.5;
        pb.scale = 0.9;
        pb.into_matrix()
    });

    chart.configure_axes().draw()?;

    chart.draw_series(
        SurfaceSeries::xoz(
            (0..=80).map(|f| f as f64 / 100.0),
            (0..=100).map(|f| f as f64 / 100.0),
            |e, ma| true_anomaly(e, ma),
        )
        .style_func(&|f| {
            let n = (f * 255.0) as u8;
            RGBColor(n, 0, 255 - n).mix(0.4).filled()
        }),
    )?;

    for e in 0..=8 {
        let e = e as f64 / 10.0;
        chart.draw_series(LineSeries::new(
            (0..=100)
                .map(|f| f as f64 / 100.0)
                .map(|ma| (e, true_anomaly(e, ma), ma)),
            &BLACK,
        ))?;
    }

    let n = 2000usize;
    let get_ma = |i: usize| i as f64 / n as f64;
    let get_ta = |i: usize, e: f64| true_anomaly(e, get_ma(i));
    let ma_iter = (0..=n).into_iter().map(get_ma);
    let ta_iter = |e: f64| ma_iter.clone().map(move |ma| true_anomaly(e, ma));

    for ta in 0..=20 {
        let target = ta as f64 / 20.0;
        chart.draw_series(LineSeries::new(
            (0..=160).map(|f| f as f64 / 200.0).map(|e| {
                let (i, ta_i) = ta_iter(e)
                    .enumerate()
                    .skip(1)
                    .find(|(_, ta)| *ta >= target)
                    .unwrap();

                let ta_i0 = get_ta(i - 1, e);
                let ma_i = get_ma(i);
                let ma_i0 = get_ma(i - 1);

                let ma_target = (target - ta_i0) / (ta_i - ta_i0) * (ma_i - ma_i0) + ma_i0;
                (e, target, ma_target)
            }),
            &BLACK,
        ))?;
    }

    chart
        .configure_series_labels()
        .border_style(&BLACK)
        .draw()?;

    // To avoid the IO failure being ignored silently, we manually call the present function
    area.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Result has been saved to {}", OUT_FILE_NAME);

    Ok(())
}

fn true_anomaly(e: f64, ma: f64) -> f64 {
    if ma < 1.0 {
        let e = Eccentricity::new(e);
        let ma = MeanAnomaly(Angle::TAU * ma);
        let ea = EccentricAnomaly::calculate(ma, e);
        let ta = TrueAnomaly::calculate(ea, e);
        ta.0 / Angle::TAU
    } else {
        1.0
    }
}
