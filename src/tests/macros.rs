macro_rules! trace_capture {
    () => {
        use tracing_subscriber::filter::{LevelFilter, Targets};
        use tracing_subscriber::prelude::*;

        let targets = Targets::new()
            .with_target("lute::tests", LevelFilter::TRACE)
            .with_target("lute::parse", LevelFilter::INFO)
            .with_target("lute::arena::molecule", LevelFilter::DEBUG)
            .with_target("lute::arena::arena", LevelFilter::TRACE);

        let formatter = tracing_subscriber::fmt::layer().with_test_writer();

        let _guard = tracing_subscriber::registry()
            .with(targets)
            .with(formatter)
            .set_default();
    };
}

pub(super) use trace_capture;
