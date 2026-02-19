use colored::Colorize;
use std::fmt::Display;
use std::time::Instant;

use crate::hardware_info::*;
use crate::stats::stat;

pub trait Measurement {
    fn new() -> Self;
    fn start(&mut self);
    fn tick(&mut self);
    fn show<T: Display>(&self, label: &str);
}

pub struct BaseTime {
    samples: Vec<f64>,
    start: Option<Instant>,
}

impl Measurement for BaseTime {
    fn new() -> Self {
        Self {
            samples: Vec::new(),
            start: None,
        }
    }

    fn start(&mut self) {
        self.samples.clear();
        self.start = Some(Instant::now());
    }

    fn tick(&mut self) {
        let start = self.start.expect("measurement not started");
        let elapsed = start.elapsed().as_secs_f64();
        self.samples.push(elapsed);
        self.start = Some(Instant::now());
    }

    fn show<T: Display>(&self, label: &str) {
        let stats = stat(&self.samples).expect("benchmark produced no samples");

        let bytes = size as f64;
        let gb = bytes / 1e9;

        let throughput_mean = gb / stats.mean;
        let throughput_stdev = gb * stats.stdev / (stats.mean * stats.mean);

        let unstable = stats.cv >= 0.10;

        let mean_str = format!("{:6.2}", throughput_mean).bright_green();

        let stdev_str = if unstable {
            format!("± {:5.2}", throughput_stdev).red().bold()
        } else {
            format!("± {:5.2}", throughput_stdev).dimmed()
        };

        let unstable_str = if unstable {
            // CV in percent
            let cv_pct = stats.cv * 100.0;
            format!(" UNSTABLE ({:.1}%)", cv_pct).red().bold()
        } else {
            "".normal()
        };

        print!("{label}:\t{} {} GB/s{}", mean_str, stdev_str, unstable_str);
        if let Some(r) = result {
            print!(" (result: {r})");
        };
        println!();
    }

    fn show_csv<T: Display>(&self, label: &str, size: u64, result: Option<T>) {
        let label = label.split_whitespace().collect::<Vec<&str>>().join(" ");
        let stats = stat(&self.samples).expect("benchmark produced no samples");

        let bytes = size as f64;
        let gb = bytes / 1e9;

        let throughput_mean = gb / stats.mean;
        let throughput_stdev = gb * stats.stdev / (stats.mean * stats.mean);

        let cpuinfo = get_hardware_info();

        print!(
            "{label},{:.6},{:.6},{:.3},{},{},{}",
            throughput_mean,
            throughput_stdev,
            stats.cv,
            cpuinfo.brand,
            cpuinfo.vendor_id,
            cpuinfo.vector_tech
        );
        if let Some(r) = result {
            print!(",{r}");
        };
        println!();
    }

    fn show_csv_header(&self, show_result: bool) {
        print!("label,mean,stdev,cv,cpu_brand,cpu_vendor,vector_ISA");
        if show_result {
            print!(",result");
        }
        println!();
    }
}
