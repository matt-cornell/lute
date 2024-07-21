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
    #[arg(short, long)]
    bench: Option<PathBuf>,
    input: String,
}
fn main() {
    let cli = Cli::parse();
    let _guard = init_tracing(cli.bench.as_deref());
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
                    OutputType::Common => match iupac_name(mol, IupacConfig::COMMON_NAME) {
                        Ok(n) => write_output(cli.out.as_deref(), n + "\n"),
                        Err(e) => tracing::error!("{e:?}"),
                    },
                    OutputType::Iupac => match iupac_name(mol, IupacConfig::PREFERRED) {
                        Ok(n) => write_output(cli.out.as_deref(), n + "\n"),
                        Err(e) => tracing::error!("{e:?}"),
                    },
                    OutputType::FastSmiles => write_output(
                        cli.out.as_deref(),
                        generate_smiles(mol, SmilesConfig::fast_roundtrip()) + "\n",
                    ),
                    OutputType::CanonSmiles => write_output(
                        cli.out.as_deref(),
                        generate_smiles(mol, SmilesConfig::new()) + "\n",
                    ),
                    OutputType::Dot => write_output(cli.out.as_deref(), fmt_as_dot(mol)),
                    OutputType::DDot => write_output(cli.out.as_deref(), AsDisp(Dot::new(mol))),
                    #[cfg(feature = "coordgen")]
                    OutputType::Svg => write_output(cli.out.as_deref(), SvgFormatter::new(mol)),
                    #[cfg(all(feature = "coordgen", feature = "resvg"))]
                    OutputType::Png => match SvgFormatter::new(mol).render(None).encode_png() {
                        Ok(b) => write_bytes(cli.out.as_deref(), &b),
                        Err(err) => tracing::error!("{err}"),
                    },
                }
            } else {
                match cli.fmt {
                    OutputType::None => {}
                    OutputType::Common => match iupac_name(&graph, IupacConfig::COMMON_NAME) {
                        Ok(n) => write_output(cli.out.as_deref(), n + "\n"),
                        Err(e) => tracing::error!("{e:?}"),
                    },
                    OutputType::Iupac => match iupac_name(&graph, IupacConfig::PREFERRED) {
                        Ok(n) => write_output(cli.out.as_deref(), n + "\n"),
                        Err(e) => tracing::error!("{e:?}"),
                    },
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
        }
        Err(err) => tracing::error!("{err}"),
    }
}
