use clap::Parser;
use lute::prelude::*;
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
    bench: Option<PathBuf>,
    #[arg(short, long)]
    rerun: Option<usize>,
    inputs: Vec<String>,
}

fn main() {
    let cli = Cli::parse();
    let _guard = init_tracing(cli.bench.as_deref());
    let mut arena = Arena::<u16>::new();
    for _ in 0..cli.rerun.unwrap_or(0) {
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
    }

    if cli.graph {
        eprintln!("arena graph: {:#?}", arena.graph());
    }
    if cli.arena {
        eprintln!("arena: {:#?}", arena.expose_frags());
    }

    let graph = lute::graph::GraphCompactor::<&UnGraph<Atom, Bond, u16>>::new(&arena.graph().inner);

    match cli.fmt {
        OutputType::None => {}
        OutputType::FastSmiles => write_output(
            cli.out.as_deref(),
            generate_smiles(&graph, SmilesConfig::fast_roundtrip()) + "\n",
        ),
        OutputType::CanonSmiles => write_output(
            cli.out.as_deref(),
            generate_smiles(&graph, SmilesConfig::new()) + "\n",
        ),
        OutputType::Dot => write_output(cli.out.as_deref(), fmt_as_dot(&graph)),
        OutputType::DDot => write_output(cli.out.as_deref(), AsDisp(Dot::new(&graph))),
        #[cfg(feature = "coordgen")]
        OutputType::Svg => write_output(cli.out.as_deref(), SvgFormatter::new(&graph)),
        #[cfg(all(feature = "coordgen", feature = "resvg"))]
        OutputType::Png => match SvgFormatter::new(&graph).render(None).encode_png() {
            Ok(b) => write_bytes(cli.out.as_deref(), &b),
            Err(err) => tracing::error!("{err}"),
        },
    }
}
