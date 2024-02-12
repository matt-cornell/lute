use clap::{Parser, ValueEnum};

#[derive(Debug, Default, Clone, Copy, ValueEnum)]
enum OutputType {
    #[default]
    Dot,
    Svg,
    Png,
}

#[derive(Parser)]
#[command(version, about)]
struct Cli {
    #[arg(short, long)]
    out: OutputType,
    input: String,
}

fn main() {
    let cli = Cli::parse();
    let parser = chem_sim::molecule::SmilesParser::new(&cli.input);
    match parser.parse() {
        Ok(graph) => match cli.out {
            OutputType::Dot => println!("{}", chem_sim::disp::fmt_as_dot(&graph)),
            OutputType::Svg => println!("{}", chem_sim::disp::fmt_as_svg(&graph)),
            OutputType::Png => todo!(),
        },
        Err(err) => eprintln!("{err}"),
    }
}
