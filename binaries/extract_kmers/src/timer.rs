use std::time::Instant;
use tracing::info;

/// Converts a `Duration` to a human-readable string like "1h23m45s"
fn human_duration(d: std::time::Duration) -> String {
    let secs = d.as_secs();
    let hours = secs / 3600;
    let minutes = (secs % 3600) / 60;
    let seconds = secs % 60;
    if hours > 0 {
        format!("{}h{:02}m{:02}s", hours, minutes, seconds)
    } else if minutes > 0 {
        format!("{}m{:02}s", minutes, seconds)
    } else {
        format!("{}s", seconds)
    }
}

/// RAII-style timer for `tracing`
pub struct TraceTimer {
    name: String,
    start: Instant,
    skippable: usize,
}

impl TraceTimer {
    pub fn new(name: impl AsRef<str>) -> Self {
        Self {
            name: name.as_ref().to_string(),
            start: Instant::now(),
            skippable: 0,
        }
    }

    pub fn skippable(name: impl AsRef<str>, skip: usize) -> Self {
        Self {
            name: name.as_ref().to_string(),
            start: Instant::now(),
            skippable: skip,
        }
    }
}

impl Drop for TraceTimer {
    fn drop(&mut self) {
        let elapsed = self.start.elapsed();
        let readable = human_duration(elapsed);
        if self.skippable <= elapsed.as_millis().try_into().unwrap() {
            info!(name = self.name, elapsed_ms = %elapsed.as_millis(), elapsed_human = %readable, "timed scope");
        }
    }
}
