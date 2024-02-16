use clap::{Parser, ValueEnum};
use std::fmt::{self, Display, Formatter};
use std::path::{Path, PathBuf};
use std::io::{Write, stdout};

#[derive(Debug, Default, Clone, Copy, ValueEnum)]
enum OutputType {
    #[default]
    Dot,
    Svg,
    Png,
}
impl Display for OutputType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dot => f.write_str("dot"),
            Self::Svg => f.write_str("svg"),
            Self::Png => f.write_str("png"),
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

fn main() {
    let cli = Cli::parse();
    let mut parser = chem_sim::parse::SmilesParser::new_unsuppressed(&cli.input);
    parser.suppress = cli.suppress;
    match parser.parse() {
        Ok(graph) => match cli.fmt {
            OutputType::Dot => write_output(cli.out.as_deref(), chem_sim::disp::fmt_as_dot(&graph)),
            OutputType::Svg => write_output(cli.out.as_deref(), chem_sim::disp::fmt_as_svg(&graph)),
            OutputType::Png => todo!(),
        },
        Err(err) => eprintln!("{err}"),
    }
}