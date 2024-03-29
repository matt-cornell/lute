use chem_sim::prelude::*;
use clap::Parser;
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
    input: String,
}
fn main() {
    init_tracing();
    let cli = Cli::parse();
    let mut parser = chem_sim::parse::SmilesParser::new(&cli.input);
    parser.suppress = !cli.unsuppress;
    match parser.parse() {
        Ok(graph) => match cli.fmt {
            OutputType::None => {}
            OutputType::Dot => write_output(cli.out.as_deref(), fmt_as_dot(&graph)),
            #[cfg(feature = "mol-svg")]
            OutputType::Svg => write_output(cli.out.as_deref(), fmt_as_svg(&graph)),
        },
        Err(err) => tracing::error!("{err}"),
    }
}
