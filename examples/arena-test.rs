use chem_sim::prelude::*;
use clap::{Parser, ValueEnum};
#[cfg(feature = "mol-bmp")]
use image::*;
use petgraph::prelude::*;
use std::fmt::{self, Display, Formatter};
#[cfg(feature = "mol-bmp")]
use std::io::{self, Cursor};
use std::io::{stdout, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutputType {
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

#[derive(Parser)]
#[command(version, about)]
struct Cli {
    #[arg(short, long, default_value_t = OutputType::None)]
    fmt: OutputType,
    #[arg(short, long)]
    out: Option<PathBuf>,
    #[arg(short, long)]
    arena: bool,
    #[arg(short, long)]
    graph: bool,
    #[arg(short, long)]
    log: bool,
    inputs: Vec<String>,
}

fn write_output<O: Display>(path: Option<&Path>, out: O) {
    let res = if let Some(p) = path {
        std::fs::File::create(p).and_then(|mut f| write!(f, "{}", out))
    } else {
        write!(stdout(), "{}", out)
    };
    if let Err(err) = res {
        eprintln!("{err}");
    }
}

fn main() {
    tracing_subscriber::fmt::init();
    let cli = Cli::parse();
    let mut arena = Arena::<u16>::new();

    for (n, input) in cli.inputs.iter().enumerate() {
        let res = SmilesParser::new(&input).parse();
        tracing::info!("parsed input {n}");
        match res {
            Ok(graph) => {
                arena.insert_mol(&graph);
                tracing::info!("inserted input {n}");
            }
            Err(err) => tracing::error!("{err}"),
        }
    }

    if cli.graph {
        eprintln!("arena graph: {:#?}", arena.graph());
    }
    if cli.arena {
        eprintln!("arena: {:#?}", arena.expose_parts());
    }

    let graph = chem_sim::graph::compact::GraphCompactor::<&StableUnGraph<Atom, Bond, u16>>::new(
        arena.graph(),
    );

    match cli.fmt {
        OutputType::None => {}
        OutputType::Dot => write_output(cli.out.as_deref(), fmt_as_dot(&graph)),
        #[cfg(feature = "mol-svg")]
        OutputType::Svg => write_output(cli.out.as_deref(), fmt_as_svg(&graph)),
    }
}
