pub struct Stat {
    pub mean: f64,
    pub stdev: f64,
    pub cv: f64,
}

pub fn stat(values: &[f64]) -> Option<Stat> {
    let n = values.len();
    if n == 0 {
        return None;
    }

    let mean = values.iter().sum::<f64>() / n as f64;

    let variance = values
        .iter()
        .map(|x| {
            let diff = x - mean;
            diff * diff
        })
        .sum::<f64>()
        / n as f64;

    let stdev = variance.sqrt();
    let cv = stdev / mean;

    Some(Stat { mean, stdev, cv })
}
