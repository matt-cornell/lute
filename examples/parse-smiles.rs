use clap::{Parser, ValueEnum};
use std::fmt::{self, Display, Formatter};
use std::io::{stdout, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Default, Clone, Copy, ValueEnum)]
enum OutputType {
    #[default]
    Show,
    Dot,
    Svg,
    Png,
}
impl Display for OutputType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Show => f.write_str("show"),
            Self::Dot => f.write_str("dot"),
            Self::Svg => f.write_str("svg"),
            Self::Png => f.write_str("png"),
        }
    }
}

#[derive(Parser)]
#[command(version, about)]
struct Cli {
    #[arg(short, long, default_value_t = OutputType::Show)]
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

fn write_png(path: Option<&Path>, img: fimg::Image<Vec<u8>, 4>) {
    if let Some(p) = path {
        match std::fs::File::create(p) {
            Ok(f) => save_png(f, img.width(), img.bytes()),
            Err(err) => eprintln!("{err}"),
        }
    } else {
        save_png(stdout(), img.width(), img.bytes());
    };
}

fn save_png<W: Write>(w: W, width: u32, data: &[u8]) {
    let mut encoder = png::Encoder::new(std::io::BufWriter::new(w), width, width);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = match encoder.write_header() {
        Ok(w) => w,
        Err(err) => {
            eprintln!("{err}");
            return
        }
    };
    if let Err(err) = writer.write_image_data(data) {
        eprintln!("{err}");
    }
}

fn main() {
    let cli = Cli::parse();
    let mut parser = chem_sim::parse::SmilesParser::new_unsuppressed(&cli.input);
    parser.suppress = cli.suppress;
    match parser.parse() {
        Ok(graph) => match cli.fmt {
        OutputType::Show => {
            if cli.out.is_some() {
                eprintln!("format was set to show, but an output file was set!");
            }
            make_img_vec(&graph).show();
        }
        OutputType::Dot => write_output(cli.out.as_deref(), fmt_as_dot(&graph)),
        OutputType::Svg => write_output(cli.out.as_deref(), fmt_as_svg(&graph)),
        OutputType::Png => write_png(cli.out.as_deref(), make_img_vec(&graph)),
    },
        Err(err) => eprintln!("{err}"),
    }
}
