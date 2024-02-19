use chem_sim::prelude::*;
use clap::{Parser, ValueEnum};
use image::*;
use std::fmt::{self, Display, Formatter};
use std::io::{self, stdout, Cursor, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutputType {
    Dot,
    Svg,
    Png,
    Jpg,
}
impl Display for OutputType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dot => f.write_str("dot"),
            Self::Svg => f.write_str("svg"),
            Self::Png => f.write_str("png"),
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
    suppress: bool,
    input: String,
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
    let mut parser = chem_sim::parse::SmilesParser::new_unsuppressed(&cli.input);
    parser.suppress = cli.suppress;
    match parser.parse() {
        Ok(graph) => match cli.fmt {
            OutputType::Dot => write_output(cli.out.as_deref(), fmt_as_dot(&graph)),
            OutputType::Svg => write_output(cli.out.as_deref(), fmt_as_svg(&graph)),
            OutputType::Png => write_out(
                cli.out.as_deref(),
                save_img(make_img_vec(&graph), ImageOutputFormat::Png),
            ),
            OutputType::Jpg => write_out(
                cli.out.as_deref(),
                save_img(make_img_vec(&graph), ImageOutputFormat::Jpeg(80)),
            ),
        },
        Err(err) => eprintln!("{err}"),
    }
}
