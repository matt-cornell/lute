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
    Dot,
    #[cfg(feature = "mol-svg")]
    Svg,
    #[cfg(feature = "mol-bmp")]
    Png,
    #[cfg(feature = "mol-bmp")]
    Jpg,
}
impl Display for OutputType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dot => f.write_str("dot"),
            #[cfg(feature = "mol-svg")]
            Self::Svg => f.write_str("svg"),
            #[cfg(feature = "mol-bmp")]
            Self::Png => f.write_str("png"),
            #[cfg(feature = "mol-bmp")]
            Self::Jpg => f.write_str("jpg"),
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

#[cfg(feature = "mol-bmp")]
fn write_out<F: FnOnce(&mut dyn Write) -> io::Result<()>>(path: Option<&Path>, callback: F) {
    let res = if let Some(p) = path {
        std::fs::File::create(p).and_then(|mut f| callback(&mut f))
    } else {
        callback(&mut stdout())
    };
    if let Err(err) = res {
        eprintln!("{err}");
    }
}

#[cfg(feature = "mol-bmp")]
fn save_img(
    img: ImageBuffer<Rgba<u8>, Vec<u8>>,
    fmt: ImageOutputFormat,
) -> impl FnOnce(&mut dyn Write) -> io::Result<()> {
    move |w| {
        let mut buf = Vec::new();
        if let Err(err) = img.write_to(&mut Cursor::new(&mut buf), fmt) {
            eprintln!("{err}");
            return Ok(());
        }
        w.write_all(&buf).map(|_| ())
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
        eprintln!("arena: {arena:#?}");
    }

    let graph = GraphCompactor::<&StableUnGraph<Atom, Bond, u16>>::new(arena.graph());

    match cli.fmt {
        OutputType::Dot => write_output(cli.out.as_deref(), fmt_as_dot(&graph)),
        #[cfg(feature = "mol-svg")]
        OutputType::Svg => write_output(cli.out.as_deref(), fmt_as_svg(&graph)),
        #[cfg(feature = "mol-bmp")]
        OutputType::Png => write_out(
            cli.out.as_deref(),
            save_img(make_img_vec(&graph), ImageOutputFormat::Png),
        ),
        #[cfg(feature = "mol-bmp")]
        OutputType::Jpg => write_out(
            cli.out.as_deref(),
            save_img(make_img_vec(&graph), ImageOutputFormat::Jpeg(80)),
        ),
    }
}
