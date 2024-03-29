use lute::prelude::*;
use clap::Parser;
use petgraph::prelude::*;
use std::path::PathBuf;

mod common;
use common::*;

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

fn main() {
    init_tracing();
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

    let graph = lute::graph::GraphCompactor::<&StableUnGraph<Atom, Bond, u16>>::new(
        arena.graph(),
    );

    match cli.fmt {
        OutputType::None => {}
        OutputType::Dot => write_output(cli.out.as_deref(), fmt_as_dot(&graph)),
        #[cfg(feature = "mol-svg")]
        OutputType::Svg => write_output(cli.out.as_deref(), fmt_as_svg(&graph)),
    }
}
