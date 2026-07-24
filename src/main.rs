//! Thin binary entry point.
//!
//! All logic lives in the `premise` library crate (`src/lib.rs`) so that the
//! alignment core can be linked by Criterion benchmarks and integration tests.
fn main() -> anyhow::Result<()> {
    premise::run()
}
