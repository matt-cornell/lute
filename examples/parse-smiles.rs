use clap::Parser;
use lute::prelude::*;
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
    unsuppress: bool,
    #[arg(short, long)]
    arena: bool,
    input: String,
}
fn main() {
    init_tracing();
    let cli = Cli::parse();
    let mut parser = SmilesParser::new(&cli.input);
    parser.suppress = !cli.unsuppress;
    match parser.parse() {
        Ok(graph) => {
            if cli.arena {
                let mut arena = Arena::<u32>::new();
                let idx = arena.insert_mol(&graph);
                let mol = arena.molecule(idx);
                match cli.fmt {
                    OutputType::None => {}
                    OutputType::Dot => write_output(cli.out.as_deref(), fmt_as_dot(mol)),
                    OutputType::DDot => write_output(cli.out.as_deref(), AsDisp(Dot::new(mol))),
                    #[cfg(feature = "mol-svg")]
                    OutputType::Svg => write_output(cli.out.as_deref(), fmt_as_svg(mol)),
                }
            } else {
                match cli.fmt {
                    OutputType::None => {}
                    OutputType::Dot => write_output(cli.out.as_deref(), fmt_as_dot(&graph)),
                    OutputType::DDot => write_output(cli.out.as_deref(), AsDisp(Dot::new(&graph))),
                    #[cfg(feature = "mol-svg")]
                    OutputType::Svg => write_output(cli.out.as_deref(), fmt_as_svg(&graph)),
                }
            }
        }
        Err(err) => tracing::error!("{err}"),
    }
}
