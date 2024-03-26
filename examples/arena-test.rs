use lute::prelude::*;
use clap::{Parser, ValueEnum};
use petgraph::prelude::*;
use std::fmt::{self, Display, Formatter};
use std::io::{stdout, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutputType {
    Dot,
    #[cfg(feature = "mol-svg")]
    Svg,
}
impl Display for OutputType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dot => f.write_str("dot"),
            #[cfg(feature = "mol-svg")]
            Self::Svg => f.write_str("svg"),
        }
    }
}

#[derive(Parser)]
#[command(version, about)]
struct Cli {
    #[arg(short, long, default_value_t = OutputType::Dot)]
    fmt: OutputType,
    #[arg(short, long)]
    out: Option<PathBuf>,
    #[arg(short, long)]
    arena: bool,
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
    let cli = Cli::parse();
    let mut arena = Arena::<u16>::new();

    for input in &cli.inputs {
        match SmilesParser::new(&input).parse() {
            Ok(graph) => {
                arena.insert_mol(&graph);
            }
            Err(err) => eprintln!("{err}"),
        }
    }

    if cli.arena {
        eprintln!("arena: {:#?}", arena.expose_parts());
    }

    let graph = GraphCompactor::<&StableUnGraph<Atom, Bond, u16>>::new(arena.graph());

    match cli.fmt {
        OutputType::Dot => write_output(cli.out.as_deref(), fmt_as_dot(&graph)),
        #[cfg(feature = "mol-svg")]
        OutputType::Svg => write_output(cli.out.as_deref(), fmt_as_svg(&graph)),
    }
}
