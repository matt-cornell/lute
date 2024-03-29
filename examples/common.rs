use clap::ValueEnum;
use std::fmt::{self, Display, Formatter};
use std::io::{stderr, stdout, Write};
use std::path::Path;

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum OutputType {
    None,
    Dot,
    #[cfg(feature = "mol-svg")]
    Svg,
}
impl Display for OutputType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::None => f.write_str("none"),
            Self::Dot => f.write_str("dot"),
            #[cfg(feature = "mol-svg")]
            Self::Svg => f.write_str("svg"),
        }
    }
}

pub fn write_output<O: Display>(path: Option<&Path>, out: O) {
    let res = if let Some(p) = path {
        std::fs::File::create(p).and_then(|mut f| write!(f, "{}", out))
    } else {
        write!(stdout(), "{}", out)
    };
    if let Err(err) = res {
        tracing::error!("{err}");
    }
}

pub fn init_tracing() {
    use tracing_subscriber::filter::*;
    use tracing_subscriber::prelude::*;
    let targets = match std::env::var("RUST_LOG") {
        Ok(var) => var.parse::<Targets>().unwrap_or_else(|e| {
            eprintln!("Ignoring `RUST_LOG={var:?}`: {e}");
            Targets::new().with_default(tracing::Level::ERROR)
        }),
        Err(e) => {
            if e != std::env::VarError::NotPresent {
                eprintln!("Ignoring `RUST_LOG`: {e}");
            }
            Targets::new().with_default(tracing::Level::ERROR)
        }
    };
    let fmt = tracing_subscriber::fmt::layer().with_writer(stderr);
    tracing_subscriber::registry()
        .with(targets)
        .with(fmt)
        .init();
}
