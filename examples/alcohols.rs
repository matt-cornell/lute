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
        OutputType::Dot => write_output(cli.out.as_deref(), fmt_as_dot(mol)),
        OutputType::DDot => write_output(cli.out.as_deref(), AsDisp(Dot::new(mol))),
        #[cfg(feature = "coordgen")]
        OutputType::CgSvg => write_output(cli.out.as_deref(), fmt_with_cg(mol)),
        #[cfg(all(feature = "coordgen", feature = "resvg"))]
        OutputType::CgPng => match fmt_with_cg(mol).render(None).encode_png() {
            Ok(b) => write_bytes(cli.out.as_deref(), &b),
            Err(err) => tracing::error!("{err}"),
        },
    }
}
