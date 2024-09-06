use clap::builder::PossibleValue;
use lute::graph::*;
use lute::prelude::*;
use petgraph::prelude::*;
use reedline_repl_rs::{self as rlr, Repl};
use rlr::clap::{self, Arg, ArgAction, ArgMatches, Command};
use std::env;
use std::ffi::OsStr;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;
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
            return &[Self::FastSmiles, Self::CanonSmiles, Self::Svg];
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
    #[error("Attempted to generate a PNG of a molecule without an output path")]
    PngToStdout,
    /// other stuff I don't wanna name
    #[error(transparent)]
    Other(#[from] Box<dyn std::error::Error + Send + Sync>),
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum FlagKind {
    Enable,
    Disable,
    Query,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum FlagName {
    Timings,
    LegacyRender,
}

#[derive(Debug, Clone, Copy)]
struct SetParser;
impl clap::builder::TypedValueParser for SetParser {
    type Value = (FlagKind, FlagName);

    fn parse_ref(
        &self,
        cmd: &Command,
        _arg: Option<&Arg>,
        value: &OsStr,
    ) -> Result<Self::Value, clap::Error> {
        use clap::error::*;
        let val = value
            .to_str()
            .ok_or_else(|| Error::new(ErrorKind::InvalidUtf8).with_cmd(cmd))?;
        let mode = match val.as_bytes().get(0) {
            Some(&b'+') => FlagKind::Enable,
            Some(&b'-') => FlagKind::Disable,
            Some(&b'?') => FlagKind::Query,
            _ => {
                let mut err = Error::new(ErrorKind::ValueValidation).with_cmd(cmd);
                err.insert(ContextKind::InvalidValue, ContextValue::String(val.to_string()));
                err.insert(
                    ContextKind::Usage,
                    ContextValue::StyledStr(
                        "Make sure the argument starts with '+', '-', or '?'".into(),
                    ),
                );
                return Err(err);
            }
        };
        let name = match &val[1..] {
            "timings" => FlagName::Timings,
            "legacy-render" | "legacy_render" => FlagName::LegacyRender,
            _ => {
                let mut err = Error::new(ErrorKind::ValueValidation).with_cmd(cmd);
                err.insert(ContextKind::InvalidValue, ContextValue::String(val.to_string()));
                err.insert(
                    ContextKind::Usage,
                    ContextValue::StyledStr(
                        r#"Available flags are "timings" and "legacy-render""#.into(),
                    ),
                );
                return Err(err);
            }
        };
        Ok((mode, name))
    }
    fn possible_values(&self) -> Option<Box<dyn Iterator<Item = PossibleValue>>> {
        Some(
            Box::new([
                PossibleValue::new("+timings").hide(true),
                PossibleValue::new("-timings").hide(true),
                PossibleValue::new("?timings"),
                PossibleValue::new("+legacy-render").alias("+legacy_render").hide(true),
                PossibleValue::new("-legacy-render").alias("-legacy_render").hide(true),
                PossibleValue::new("?legacy-render").alias("?legacy_render"),
            ]
            .into_iter()),
        )
    }
}

#[derive(Debug, Clone)]
struct Timer(Instant);
impl Timer {
    pub fn new() -> Self {
        Self(Instant::now())
    }
}
impl Drop for Timer {
    fn drop(&mut self) {
        println!("{:?}", self.0.elapsed());
    }
}

type SubscriberType = tracing_subscriber::layer::Layered<
    Targets,
    fmt::Subscriber<DefaultFields, Format<Full>, LevelFilter, fn() -> std::io::Stderr>,
>;

struct Context {
    arena: Arena<u32>,
    trace: Arc<SubscriberType>,
    timings: bool,
    legacy_render: bool,
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
    if let Some(filter) = args.get_one::<String>("filter") {
        ctx.trace = Arc::new(make_subscriber(filter)?);
        Ok(None)
    } else {
        Ok(Some(
            ctx.trace.downcast_ref::<Targets>().unwrap().to_string(),
        ))
    }
}

fn load_smiles(args: ArgMatches, ctx: &mut Context) -> Result<Option<String>, ReplError> {
    use std::fmt::Write;
    let _guard = tracing::subscriber::set_default(ctx.trace.clone());
    let _timer = ctx.timings.then(Timer::new);
    let mut out = String::new();
    let inputs = args.get_many::<String>("input").unwrap();
    let len = inputs.len();
    for (n, arg) in inputs.enumerate() {
        let graph = SmilesParser::new(arg)
            .parse()
            .map_err(|err| ReplError::Smiles {
                err,
                index: if len == 1 { 0 } else { n },
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
    let _guard = tracing::subscriber::set_default(ctx.trace.clone());
    let _timer = ctx.timings.then(Timer::new);
    let path = args.get_one::<PathBuf>("out");
    match args.get_one::<DumpType>("type").unwrap() {
        ty @ (DumpType::CanonSmiles | DumpType::FastSmiles) => {
            let cfg = if *ty == DumpType::FastSmiles {
                SmilesConfig::fast_roundtrip()
            } else {
                SmilesConfig::new()
            };
            let s = if let Some(inputs) = args.get_many::<u32>("input") {
                generate_smiles(
                    &GraphUnion(inputs.map(|i| ctx.arena.molecule(*i)).collect()),
                    cfg,
                )
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
        #[cfg(feature = "coordgen")]
        DumpType::Svg => {
            use lute::disp::svg::FormatMode;
            let mode = if ctx.legacy_render {
                FormatMode::LegacyH
            } else {
                FormatMode::Normal
            };
            let s = if let Some(inputs) = args.get_many::<u32>("input") {
                SvgFormatter {
                    graph: &GraphUnion(inputs.map(|i| ctx.arena.molecule(*i)).collect()),
                    mode,
                }
                .to_string()
            } else {
                SvgFormatter {
                    graph: &GraphCompactor::<&StableUnGraph<Atom, Bond>>::new(ctx.arena.graph()),
                    mode,
                }
                .to_string()
            };
            if let Some(p) = path {
                std::fs::write(&p, s)?;
                Ok(Some(format!("saved to {}", p.display())))
            } else {
                Ok(Some(s))
            }
        }
        #[cfg(all(feature = "coordgen", feature = "resvg"))]
        DumpType::Png => {
            use lute::disp::svg::FormatMode;
            let mode = if ctx.legacy_render {
                FormatMode::LegacyH
            } else {
                FormatMode::Normal
            };
            let s = if let Some(inputs) = args.get_many::<u32>("input") {
                SvgFormatter {
                    graph: &GraphUnion(inputs.map(|i| ctx.arena.molecule(*i)).collect()),
                    mode,
                }
                .render(None)
            } else {
                SvgFormatter {
                    graph: &GraphCompactor::<&StableUnGraph<Atom, Bond>>::new(ctx.arena.graph()),
                    mode,
                }
                .render(None)
            };
            if let Some(p) = path {
                s.save_png(&p).map_err(Box::from)?;
                Ok(Some(format!("saved to {}", p.display())))
            } else {
                Err(ReplError::PngToStdout)
            }
        }
    }
}

fn adjust_settings(args: ArgMatches, ctx: &mut Context) -> Result<Option<String>, ReplError> {
    use std::fmt::Write;
    let mut out = String::new();
    for (act, name) in args.get_many::<(FlagKind, FlagName)>("settings").unwrap() {
        let flag = match name {
            FlagName::Timings => &mut ctx.timings,
            FlagName::LegacyRender => &mut ctx.legacy_render,
        };
        match act {
            FlagKind::Enable => *flag = true,
            FlagKind::Disable => *flag = false,
            FlagKind::Query => {
                if out.is_empty() {
                    out.push(' ');
                    let _ = write!(out, "{flag}");
                }
            }
        }
    }
    Ok(Some(out))
}

fn main() {
    let trace = match env::var("LUTE_LOG") {
        Ok(name) => make_subscriber(&name).unwrap_or_else(|e| {
            eprintln!("Error in LUTE_LOG environment variable: {e}");
            default_subscriber()
        }),
        Err(env::VarError::NotPresent) => default_subscriber(),
        Err(env::VarError::NotUnicode(_)) => {
            eprintln!("Error in LUTE_LOG environment variable: not valid UTF-8");
            default_subscriber()
        }
    };
    let context = Context {
        arena: Arena::new(),
        trace: Arc::new(trace),
        timings: false,
        legacy_render: false,
    };
    let mut repl = Repl::new(context)
        .with_name("lute")
        .with_version(concat!("v", env!("CARGO_PKG_VERSION")))
        .with_prompt(concat!("lute v", env!("CARGO_PKG_VERSION")))
        .with_command(Command::new("reset").about("Reset the arena"), |_, ctx| {
            ctx.arena = Arena::new();
            Ok(None)
        })
        .with_command(
            Command::new("log")
                .about("Set the log filter")
                .arg(Arg::new("filter").required(false).help("The filter string to use")),
            set_log,
        )
        .with_command(
            Command::new("load")
                .about("Load SMILES inputs")
                .arg(Arg::new("input").required(true).action(ArgAction::Append).help("SMILES inputs to load")),
            load_smiles,
        )
        .with_command(
            Command::new("dump")
                .about("Show the output of inputs, or the whole graph if no ids are given")
                .arg(
                    Arg::new("type")
                        .required(true)
                        .value_parser(clap::value_parser!(DumpType))
                        .help("Output format to dump as"),
                )
                .arg(
                    Arg::new("out")
                        .short('o')
                        .long("output")
                        .value_name("FILE")
                        .value_parser(clap::value_parser!(PathBuf))
                        .help("The file to output to"),
                )
                .arg(
                    Arg::new("input")
                        .required(false)
                        .action(ArgAction::Append)
                        .value_name("IDS")
                        .value_parser(clap::value_parser!(u32))
                        .help("Fragment IDs to dump"),
                ),
            dump,
        )
        .with_command(
            Command::new("set").about("Adjust REPL settings").arg(
                Arg::new("settings")
                    .required(true)
                    .action(ArgAction::Append)
                    .allow_hyphen_values(true)
                    .value_parser(clap::builder::ValueParser::new(SetParser)),
            ),
            adjust_settings,
        );
    let res = repl.run();
    if let Err(err) = res {
        eprintln!("{err}");
    }
}
