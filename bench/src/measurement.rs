use std::fmt::Display;
use std::time::Instant;

pub trait Measurement {
    fn start(&mut self);
    fn show<T: Display>(&mut self, label: &str, size: u64, rep: u64, result: T);
    fn new() -> Self;
}

pub struct BaseTime(Option<Instant>);

impl Measurement for BaseTime {
    fn new() -> Self {
        Self(None)
    }
    fn start(&mut self) {
        self.0 = Some(Instant::now());
    }
    fn show<T: Display>(&mut self, label: &str, size: u64, rep: u64, result: T) {
        let val = self.0.unwrap().elapsed().as_secs_f64();
        println!(
            "{label}:\t {:5.2} GB/s  result:{}",
            (size * rep) as f64 / 1e9 / val,
            result
        );
    }
}
