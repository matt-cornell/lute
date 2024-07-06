use clap::ValueEnum;
use std::fmt::{self, Debug, Display, Formatter};
use std::io::{stderr, stdout, Write};
use std::path::Path;

pub use petgraph::dot::Dot;

#[derive(Debug, Clone, Copy)]
pub struct AsDisp<T>(pub T);
impl<T: Debug> Display for AsDisp<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        Debug::fmt(&self.0, f)
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum OutputType {
    None,
    FastSmiles,
    CanonSmiles,
    Dot,
    DDot,
    #[cfg(feature = "mol-svg")]
    Svg,
    #[cfg(all(feature = "mol-svg", feature = "resvg"))]
    Png,
}
impl Display for OutputType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::None => f.write_str("none"),
            Self::FastSmiles => f.write_str("fast-smiles"),
            Self::CanonSmiles => f.write_str("canon-smiles"),
            Self::Dot => f.write_str("dot"),
            Self::DDot => f.write_str("d-dot"),
            #[cfg(feature = "mol-svg")]
            Self::Svg => f.write_str("svg"),
            #[cfg(all(feature = "mol-svg", feature = "resvg"))]
            Self::Png => f.write_str("png"),
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

#[allow(dead_code)]
pub fn write_bytes(path: Option<&Path>, bytes: &[u8]) {
    let res = if let Some(p) = path {
        std::fs::File::create(p).and_then(|mut f| f.write_all(bytes))
    } else {
        stdout().write_all(bytes)
    };
    if let Err(err) = res {
        tracing::error!("{err}");
    }
}

pub fn init_tracing(flame: Option<&Path>) -> Option<tracing_flame::FlushGuard<impl Write>> {
    use tracing_subscriber::filter::*;
    use tracing_subscriber::prelude::*;
    use tracing_subscriber::*;
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

    if let Some(flame) = flame {
        let fmt = fmt::layer().with_writer(stderr);
        let (flame, guard) =
            tracing_flame::FlameLayer::with_file(flame).expect("FlameLayer failed");
        registry()
            .with(fmt.with_filter(targets))
            .with(flame.with_module_path(false))
            .init();
        Some(guard)
    } else {
        let fmt = fmt::fmt().with_writer(stderr).finish();
        fmt.with(targets).init();
        None
    }
}
