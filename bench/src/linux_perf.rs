use colored::*;
use perf_event::events::Hardware;
use perf_event::{Builder, Counter};
use std::fmt::Display;
use std::time::Instant; // terminal colors

use crate::measurement::Measurement;

fn human_readable(n: u64) -> String {
    // returns 1.23K, 45.6M, 7.89B, etc.
    let n_f = n as f64;
    if n < 1_000 {
        format!("{}", n)
    } else if n < 1_000_000 {
        format!("{:.2}K", n_f / 1_000.0)
    } else if n < 1_000_000_000 {
        format!("{:.2}M", n_f / 1_000_000.0)
    } else {
        format!("{:.2}G", n_f / 1_000_000_000.0)
    }
}

pub struct PerfMeasurement {
    start: Option<Instant>,
    cycles: Counter,
    instructions: Counter,
    branches: Counter,
    branch_misses: Counter,

    // store readings after stop
    cycles_val: u64,
    instructions_val: u64,
    branches_val: u64,
    branch_misses_val: u64,
}

impl Measurement for PerfMeasurement {
    fn new() -> Self {
        let mut cycles = Builder::new()
            .kind(Hardware::CPU_CYCLES)
            .build()
            .expect("failed to create perf counter for cycles");
        cycles.disable().unwrap();

        let mut instructions = Builder::new()
            .kind(Hardware::INSTRUCTIONS)
            .build()
            .expect("failed to create perf counter for instructions");
        instructions.disable().unwrap();

        let mut branches = Builder::new()
            .kind(Hardware::BRANCH_INSTRUCTIONS)
            .build()
            .expect("failed to create perf counter for branches");
        branches.disable().unwrap();

        let mut branch_misses = Builder::new()
            .kind(Hardware::BRANCH_MISSES)
            .build()
            .expect("failed to create perf counter for branch misses");
        branch_misses.disable().unwrap();

        Self {
            start: None,
            cycles,
            instructions,
            branches,
            branch_misses,
            cycles_val: 0,
            instructions_val: 0,
            branches_val: 0,
            branch_misses_val: 0,
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
    }

    fn show<T: Display>(&mut self, label: &str, size: u64, rep: u64, result: T) {
        let elapsed = self.start.unwrap().elapsed().as_secs_f64();
        let bytes_total = (size * rep) as f64;

        // read counters
        self.cycles_val = self.cycles.read().unwrap();
        self.instructions_val = self.instructions.read().unwrap();
        self.branches_val = self.branches.read().unwrap();
        self.branch_misses_val = self.branch_misses.read().unwrap();

        let cycles_per_byte = self.cycles_val as f64 / bytes_total;
        let instr_per_byte = self.instructions_val as f64 / bytes_total;
        let branches_per_byte = self.branches_val as f64 / bytes_total;
        let branch_miss_pct = if self.branches_val > 0 {
            100.0 * self.branch_misses_val as f64 / self.branches_val as f64
        } else {
            0.0
        };
        let gbps = bytes_total / 1e9 / elapsed;

        // human-readable numbers
        let cycles_str = human_readable(self.cycles_val).bright_blue();
        let instr_str = human_readable(self.instructions_val).bright_magenta();
        let branches_str = human_readable(self.branches_val).bright_yellow();
        let branch_miss_str = human_readable(self.branch_misses_val).bright_red();

        // floats
        let gbps_str = format!("{:>8.2}GB/s", gbps).bright_green();
        let cpb_str = format!("{:>8.2}", cycles_per_byte).bright_blue();
        let ipb_str = format!("{:>8.2}", instr_per_byte).bright_magenta();
        let bpb_str = format!("{:>8.2}", branches_per_byte).bright_yellow();
        let branch_miss_pct_str = format!("{:>6.2}%", branch_miss_pct).bright_red();

        // aligned labels
        println!("");
        println!("{}:", label.bold());
        println!("    {:>16} : {}", "throughput", gbps_str);
        println!("    {:>16} : {}", "cycles", cycles_str);
        println!("    {:>16} : {}", "cycles/byte", cpb_str);
        println!("    {:>16} : {}", "instructions", instr_str);
        println!("    {:>16} : {}", "instr/byte", ipb_str);
        println!("    {:>16} : {}", "branches", branches_str);
        println!("    {:>16} : {}", "branches/byte", bpb_str);
        println!("    {:>16} : {}", "% branch miss", branch_miss_pct_str);
        println!("    {:>16} : {}", "branch misses", branch_miss_str);
        println!("    {:>16} : {}", "result", result);
    }
}
