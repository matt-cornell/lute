use lute::prelude::*;
use reedline_repl_rs::{self as rlr, Repl};
use rlr::clap::{self, Arg, ArgAction, ArgMatches, Command};
use clap::builder::PossibleValue;
use std::convert::Infallible;
use std::path::PathBuf;
use std::env;
use thiserror::Error;
use tracing_subscriber::filter::{LevelFilter, ParseError, Targets};
use tracing_subscriber::fmt::{self, format::*};
use tracing_subscriber::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq)]
enum DumpType {
    FastSmiles,
    CanonSmiles,
    #[cfg(feature = "coordgen")]
    Svg,
    #[cfg(all(feature = "coordgen", feature = "resvg"))]
    Png,
}
impl clap::ValueEnum for DumpType {
    #[allow(unreachable_code)]
    fn value_variants<'a>() -> &'a [Self] {
        #[cfg(feature = "coordgen")]
        {
            #[cfg(feature = "resvg")]
            return &[Self::FastSmiles, Self::CanonSmiles, Self::Svg, Self::Png];
            return &[Self::FastSmiles, Self::CanonSmiles, Self::Svg]
        }
        &[Self::FastSmiles, Self::CanonSmiles]
    }
    fn to_possible_value(&self) -> Option<PossibleValue> {
        match self {
            Self::FastSmiles => Some(PossibleValue::new("fast-smiles")),
            Self::CanonSmiles => Some(PossibleValue::new("canon-smiles").alias("smiles")),
            #[cfg(feature = "coordgen")]
            Self::Svg => Some(PossibleValue::new("svg")),
            #[cfg(all(feature = "coordgen", feature = "resvg"))]
            Self::Png => Some(PossibleValue::new("png")),
        }
    }

}

#[derive(Debug, Error)]
enum ReplError {
    #[error(transparent)]
    Repl(#[from] rlr::Error),
    #[error(transparent)]
    Tracing(#[from] ParseError),
    #[error(transparent)]
    StdIo(#[from] std::io::Error),
    #[error("{}{err}", if *.index == 0 { String::new() } else { format!("Error in input #{index}: ") })]
    Smiles {
        index: usize,
        err: lute::parse::smiles::SmilesError,
    },
}

type SubscriberType = tracing_subscriber::layer::Layered<
    Targets,
    fmt::Subscriber<DefaultFields, Format<Full>, LevelFilter, fn() -> std::io::Stderr>,
>;

struct Context {
    arena: Arena<u32>,
    trace: SubscriberType,
}

fn default_subscriber() -> SubscriberType {
    fmt::fmt()
        .with_writer(std::io::stderr as _)
        .finish()
        .with(Targets::new())
}

fn make_subscriber(filter: &str) -> Result<SubscriberType, ParseError> {
    use std::str::FromStr;
    Ok(fmt::fmt()
        .with_writer(std::io::stderr as _)
        .finish()
        .with(Targets::from_str(filter)?))
}

fn set_log(args: ArgMatches, ctx: &mut Context) -> Result<Option<String>, ReplError> {
    ctx.trace = make_subscriber(args.get_one::<String>("filter").unwrap())?;
    Ok(None)
}

fn load_smiles(args: ArgMatches, ctx: &mut Context) -> Result<Option<String>, ReplError> {
    use std::fmt::Write;
    let mut out = String::new();
    let inputs = args.get_many::<String>("input").unwrap();
    let len = inputs.len();
    for (n, arg) in inputs.enumerate() {
        let graph = SmilesParser::new(arg)
            .parse()
            .map_err(|err| ReplError::Smiles {
                err,
                index: if len == 1 { 0 } else { len },
            })?;
        let id = ctx.arena.insert_mol(&graph);
        if !out.is_empty() {
            out.push(' ');
        }
        let _ = write!(out, "{id}");
    }
    Ok(Some(out))
}

fn dump(args: ArgMatches, ctx: &mut Context) -> Result<Option<String>, ReplError> {
    let path = args.get_one::<PathBuf>("out");
    match args.get_one::<DumpType>("type").unwrap() {
        ty @ (DumpType::CanonSmiles | DumpType::FastSmiles) => {
            let cfg = if *ty == DumpType::FastSmiles {
                SmilesConfig::fast_roundtrip()
            } else {
                SmilesConfig::new()
            };
            let s = if let Some(inputs) = args.get_many::<u32>("input") {
                todo!()
            } else {
                generate_smiles(ctx.arena.graph(), cfg)
            };
            if let Some(p) = path {
                std::fs::write(&p, s)?;
                Ok(Some(format!("saved to {}", p.display())))
            } else {
                Ok(Some(s))
            }
        }
        _ => unreachable!(),
    }
}

fn main() {
    let trace = match env::var("LUTE_LOG") {
        Ok(name) => {
            make_subscriber(&name).unwrap_or_else(|e| {
                eprintln!("Error in LUTE_LOG environment variable: {e}");
                default_subscriber()
            })
        }
        Err(env::VarError::NotPresent) => default_subscriber(),
        Err(env::VarError::NotUnicode(_)) => {
            eprintln!("Error in LUTE_LOG environment variable: not valid UTF-8");
            default_subscriber()
        }
    };
    let context = Context {
        arena: Arena::new(),
        trace,
    };
    let mut repl = Repl::new(context)
        .with_name("Lute")
        .with_command(Command::new("reset").about("Reset the arena"), |_, ctx| {
            ctx.arena = Arena::new();
            Ok(None)
        })
        .with_command(
            Command::new("log")
                .about("Set the log filter")
                .arg(Arg::new("filter").required(true)),
            set_log,
        )
        .with_command(
            Command::new("load")
                .about("Load SMILES inputs")
                .arg(Arg::new("input").required(true).action(ArgAction::Append)),
            load_smiles,
        )
        .with_command(
            Command::new("dump")
                .about("Show the output of inputs, or the whole graph if no ids are given")
                .arg(Arg::new("type").required(true).value_parser(clap::value_parser!(DumpType)))
                .arg(Arg::new("out").short('o').long("output").value_parser(clap::value_parser!(PathBuf)))
                .arg(
                    Arg::new("input")
                        .required(false)
                        .action(ArgAction::Append)
                        .value_parser(clap::value_parser!(u32)),
                ),
            dump,
        );
    repl.run();
}
