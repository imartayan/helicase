fn main() {
    // When invoked via `cargo test --all-targets`, cargo executes benchmark binaries
    // in unoptimized debug mode. Exit immediately so test suites finish in seconds.
    if cfg!(debug_assertions)
        && !std::env::args().any(|a| a == "--bench" || a == "-r" || a == "--repeat")
    {
        return;
    }
    bench::cli::run();
}
