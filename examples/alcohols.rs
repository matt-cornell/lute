use clap::Parser;
use lute::prelude::*;
use std::path::PathBuf;

mod common;
use common::*;

#[derive(Parser)]
struct Cli {
    #[arg(short, long, default_value_t = OutputType::None)]
    fmt: OutputType,
    #[arg(short, long)]
    out: Option<PathBuf>,
    #[arg(short, long)]
    bench: Option<PathBuf>,
}

fn main() {
    let cli = Cli::parse();
    let _guard = init_tracing(cli.bench.as_deref());

    let mut arena = Arena::<u8>::new();

    let secondary_alcohol = smiles!("C&&O");
    arena.insert_mol(&secondary_alcohol);

    let isopropanol = smiles!("CC(C)O");
    let isp_idx = arena.insert_mol(&isopropanol);

    let mol = arena.molecule(isp_idx);

    match cli.fmt {
        OutputType::None => {}
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
        #[cfg(feature = "mol-svg")]
        OutputType::Svg => write_output(cli.out.as_deref(), fmt_as_svg(mol)),
        #[cfg(all(feature = "mol-svg", feature = "resvg"))]
        OutputType::Png => match fmt_as_svg(mol).render(None).encode_png() {
            Ok(b) => write_bytes(cli.out.as_deref(), &b),
            Err(err) => tracing::error!("{err}"),
        },
    }
}
