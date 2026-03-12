use colored::*;
use helicase::PDEP_ACTIVATE;
use perf_event::events::Hardware;
use perf_event::{Builder, Counter};
use std::fmt::Display;
use std::time::Instant;

use crate::hardware_info::*;
use crate::measurement::Measurement;
use crate::stats::*;

fn human_readable(n: f64) -> String {
    if n < 1_000.0 {
        format!("{:.2}", n)
    } else if n < 1_000_000.0 {
        format!("{:.2}K", n / 1_000.0)
    } else if n < 1_000_000_000.0 {
        format!("{:.2}M", n / 1_000_000.0)
    } else {
        format!("{:.2}G", n / 1_000_000_000.0)
    }
}

/// Format a metric (mean ± stdev) and highlight if unstable
fn fmt_stat(stat: &Stat, scale: f64, threshold: f64, unit: Option<&str>) -> colored::ColoredString {
    let mean = human_readable(stat.mean * scale);
    let mut s = mean;
    if let Some(u) = unit {
        s = format!("{} {}", s, u);
    }
    s = format!("{:>10} ± {}", s.bold(), human_readable(stat.stdev * scale));
    if stat.cv >= threshold {
        format!("{} UNSTABLE", s).red().bold()
    } else {
        s.normal()
    }
}

pub struct PerfMeasurement {
    filename: String,
    start: Option<Instant>,
    cycles: Counter,
    instructions: Counter,
    branches: Counter,
    branch_misses: Counter,

    // store readings after stop
    time: Vec<f64>,
    cycles_val: Vec<u64>,
    instructions_val: Vec<u64>,
    branches_val: Vec<u64>,
    branch_misses_val: Vec<u64>,
}

impl Measurement for PerfMeasurement {
    fn new(filename: &str) -> Self {
        let mut cycles = Builder::new(Hardware::CPU_CYCLES)
            .build()
            .expect("failed to create perf counter for cycles");
        cycles.disable().unwrap();

        let mut instructions = Builder::new(Hardware::INSTRUCTIONS)
            .build()
            .expect("failed to create perf counter for instructions");
        instructions.disable().unwrap();

        let mut branches = Builder::new(Hardware::BRANCH_INSTRUCTIONS)
            .build()
            .expect("failed to create perf counter for branches");
        branches.disable().unwrap();

        let mut branch_misses = Builder::new(Hardware::BRANCH_MISSES)
            .build()
            .expect("failed to create perf counter for branch misses");
        branch_misses.disable().unwrap();

        Self {
            filename: filename.to_string(),
            start: None,
            cycles,
            instructions,
            branches,
            branch_misses,
            time: vec![],
            cycles_val: vec![],
            instructions_val: vec![],
            branches_val: vec![],
            branch_misses_val: vec![],
        }
    }

    fn start(&mut self) {
        self.start = Some(Instant::now());

        self.cycles.reset().unwrap();
        self.instructions.reset().unwrap();
        self.branches.reset().unwrap();
        self.branch_misses.reset().unwrap();

        self.cycles.enable().unwrap();
        self.instructions.enable().unwrap();
        self.branches.enable().unwrap();
        self.branch_misses.enable().unwrap();

        self.time.clear();
        self.cycles_val.clear();
        self.instructions_val.clear();
        self.branches_val.clear();
        self.branch_misses_val.clear();
    }

    fn tick(&mut self) {
        let start = self.start.expect("measurement not started");
        let elapsed = start.elapsed().as_secs_f64();
        self.time.push(elapsed);
        self.start = Some(Instant::now()); // reset timer for next tick

        self.cycles_val.push(self.cycles.read().unwrap());
        self.instructions_val
            .push(self.instructions.read().unwrap());
        self.branches_val.push(self.branches.read().unwrap());
        self.branch_misses_val
            .push(self.branch_misses.read().unwrap());

        // reset counters after reading
        self.cycles.reset().unwrap();
        self.instructions.reset().unwrap();
        self.branches.reset().unwrap();
        self.branch_misses.reset().unwrap();
    }

    fn show_human<T: Display>(&self, label: &str, size: u64, result: Option<T>) {
        let bytes = size as f64;
        let gb = bytes / 1e9;
        let threshold = 0.10; // CV >= 10% considered unstable

        // Convert counters to f64 once
        let cycles_f: Vec<f64> = self.cycles_val.iter().map(|&v| v as f64).collect();
        let instr_f: Vec<f64> = self.instructions_val.iter().map(|&v| v as f64).collect();
        let branches_f: Vec<f64> = self.branches_val.iter().map(|&v| v as f64).collect();
        let branch_misses_f: Vec<f64> = self.branch_misses_val.iter().map(|&v| v as f64).collect();

        // Compute stats
        let cycles = stat(&cycles_f).expect("no cycle samples");
        let instr = stat(&instr_f).expect("no instruction samples");
        let branches = stat(&branches_f).expect("no branch samples");
        let branch_misses = stat(&branch_misses_f).expect("no branch-miss samples");
        let time = stat(&self.time).expect("no time samples");

        // Throughput GB/s
        let gbps_mean = gb / time.mean;
        let gbps_stdev = gb * time.stdev / (time.mean * time.mean);
        let gbps_str = fmt_stat(
            &Stat {
                mean: gbps_mean,
                stdev: gbps_stdev,
                cv: time.cv,
            },
            1.0,
            threshold,
            Some("GB/s"),
        );

        let cpb_str = fmt_stat(&cycles, 1.0 / bytes, threshold, None);
        let ipb_str = fmt_stat(&instr, 1.0 / bytes, threshold, None);
        let bpb_str = fmt_stat(&branches, 1.0 / bytes, threshold, None);

        let branch_miss_pct = Stat {
            mean: branch_misses.mean / branches.mean,
            stdev: 0.0,
            cv: branch_misses.cv,
        };
        let branch_miss_str = fmt_stat(&branch_miss_pct, 1000.0, threshold, Some("‰"));

        // Human-readable totals
        let cycles_str = fmt_stat(&cycles, 1.0, threshold, None).to_string();
        let instr_str = fmt_stat(&instr, 1.0, threshold, None).to_string();
        let branches_str = fmt_stat(&branches, 1.0, threshold, None).to_string();
        let branch_misses_str = fmt_stat(&branch_misses, 1.0, threshold, None).to_string();

        // Print aligned
        println!("\n{}:", label.bold());
        println!("    {:>16} : {}", "throughput", gbps_str);
        println!("    {:>16} : {}", "cycles", cycles_str);
        println!("    {:>16} : {}", "cycles/byte", cpb_str);
        println!("    {:>16} : {}", "instructions", instr_str);
        println!("    {:>16} : {}", "instr/byte", ipb_str);
        println!("    {:>16} : {}", "branches", branches_str);
        println!("    {:>16} : {}", "branches/byte", bpb_str);
        println!("    {:>16} : {}", "‰ branch miss", branch_miss_str);
        println!("    {:>16} : {}", "branch misses", branch_misses_str);
        if let Some(r) = result {
            println!("    {:>16} : {r}", "result");
        }
    }

    fn show_csv_header(&self, show_result: bool) {
        print!(
            "label,time_mean,time_stdev,cycle_mean,cycle_stdev,instructions_mean,instructions_stdev,branches_mean,branches_stdev,branch_misses_mean,branch_misses_stdev,size,filename,cpu_brand,cpu_vendor,vector_ISA,pdep_activate"
        );
        if show_result {
            print!(",result");
        }
        println!();
    }

    fn show_csv<T: Display>(&self, label: &str, size: u64, result: Option<T>) {
        let label = label.split_whitespace().collect::<Vec<&str>>().join(" ");

        // Convert counters to f64 once
        let cycles_f: Vec<f64> = self.cycles_val.iter().map(|&v| v as f64).collect();
        let instr_f: Vec<f64> = self.instructions_val.iter().map(|&v| v as f64).collect();
        let branches_f: Vec<f64> = self.branches_val.iter().map(|&v| v as f64).collect();
        let branch_misses_f: Vec<f64> = self.branch_misses_val.iter().map(|&v| v as f64).collect();

        // Compute stats
        let cycles = stat(&cycles_f).expect("no cycle samples");
        let instr = stat(&instr_f).expect("no instruction samples");
        let branches = stat(&branches_f).expect("no branch samples");
        let branch_misses = stat(&branch_misses_f).expect("no branch-miss samples");
        let time = stat(&self.time).expect("no time samples");
        let cpuinfo = get_hardware_info();

        print!(
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            label,
            time.mean,
            time.stdev,
            cycles.mean,
            cycles.stdev,
            instr.mean,
            instr.stdev,
            branches.mean,
            branches.stdev,
            branch_misses.mean,
            branch_misses.stdev,
            size,
            self.filename,
            cpuinfo.brand,
            cpuinfo.vendor_id,
            cpuinfo.vector_tech,
            PDEP_ACTIVATE,
        );
        if let Some(r) = result {
            print!(",{r}");
        };
        println!();
    }
}
