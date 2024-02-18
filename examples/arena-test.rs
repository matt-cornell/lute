use clap::{Parser, ValueEnum};
use std::fmt::{self, Display, Formatter};
use std::io::{stdout, Write};
use std::path::{Path, PathBuf};
use chem_sim::prelude::*;
use petgraph::prelude::*;

#[derive(Debug, Default, Clone, Copy, ValueEnum)]
enum OutputType {
    #[default]
    Dot,
    Svg,
    Png,
}
impl Display for OutputType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dot => f.write_str("dot"),
            Self::Svg => f.write_str("svg"),
            Self::Png => f.write_str("png"),
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

    let graph = GraphCompactor::<&StableUnGraph<Atom, Bond, u16>>::new(arena.graph());

    match cli.fmt {
        OutputType::Dot => write_output(cli.out.as_deref(), fmt_as_dot(&graph)),
        OutputType::Svg => write_output(cli.out.as_deref(), fmt_as_svg(&graph)),
        OutputType::Png => todo!(),
    }
}
