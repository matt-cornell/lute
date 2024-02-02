use chem_sim::atom_info::*;
use chem_sim::molecule::*;
use clap::{Parser, ValueEnum};
use petgraph::dot::*;
use petgraph::graph::*;
use petgraph::visit::EdgeRef;

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

fn dot_atom_color(num: u8) -> &'static str {
    match num {
        0 => "chartreuse",
        1 => "lightblue2",
        3 => "purple1",
        4 => "magenta1",
        5 => "darkorange4",
        6 => "gray16",
        7 => "deepskyblue",
        8 => "red",
        9 => "green1",
        11 => "purple2",
        12 => "magenta2",
        14 => "khaki4",
        15 => "firebrick3",
        16 => "gold3",
        17 => "green2",
        19 => "purple3",
        20 => "magenta3",
        35 => "green4",

        _ => match ATOM_DATA[num as usize].group {
            ElemGroup::Nonmet => "darkorchid4",
            ElemGroup::Noble => "lightcoral",
            ElemGroup::Trans => "slategray2",
            ElemGroup::Poor => "slategray4",
            ElemGroup::RarEar => "maroon",
            ElemGroup::Metoid => "lightgoldenrod4",
            ElemGroup::Alkali => "purple4",
            ElemGroup::AlkEar => "magenta4",
            ElemGroup::Halogn => "green4",
        },
    }
}
fn svg_atom_color(num: u8) -> &'static str {
    match num {
        0 => "#7FFF00",
        1 => "#ADD8E6",
        3 => "#9400D3",
        4 => "#8B008B",
        5 => "#CC7000",
        6 => "#050505",
        7 => "#00BFFF",
        8 => "#CC0000",
        9 => "#00C000",
        11 => "#7900AD",
        12 => "#750075",
        14 => "khaki4",
        15 => "#8B0000",
        16 => "#DAA520",
        17 => "#00A000",
        19 => "#60008A",
        20 => "#610061",
        35 => "#008000",

        _ => match ATOM_DATA[num as usize].group {
            ElemGroup::Nonmet => "#20B2AA",
            ElemGroup::Noble => "#FF69B4",
            ElemGroup::Trans => "#778899",
            ElemGroup::Poor => "#2F4F4F",
            ElemGroup::RarEar => "#FF1943",
            ElemGroup::Metoid => "#556B2F",
            ElemGroup::Alkali => "#470066",
            ElemGroup::AlkEar => "#610061",
            ElemGroup::Halogn => "#006000",
        },
    }
}
fn bond_style<'a>(_: &'a MoleculeGraph, bond: EdgeReference<'a, Bond>) -> String {
    match bond.weight() {
        Bond::Single => r#"color="black""#,
        Bond::Double => r#"color="black:black""#,
        Bond::Triple => r#"color="black:black:black""#,
        Bond::Quad => r#"color="black:black:black:black""#,
        Bond::Aromatic => r#"color="black:black", style="dashed""#,
        Bond::Left => r#"color="black", label="left""#,
        Bond::Right => r#"color="black", label="right""#,
        Bond::Non => r#"color="black", style="invis""#,
    }
    .to_string()
}
fn atom_style<'a>(_: &'a MoleculeGraph, (_, atom): (NodeIndex, &'a Atom)) -> String {
    let num = atom.protons;
    format!("color={:?}", dot_atom_color(num))
}

fn main() {
    let cli = Cli::parse();
    let parser = SmilesParser::new(&cli.input);
    match parser.parse() {
        Ok(graph) => match cli.out {
            OutputType::Dot => println!(
                "{}",
                Dot::with_attr_getters(&graph, &[Config::EdgeNoLabel], &bond_style, &atom_style)
            ),
            OutputType::Svg => {
                let atoms = graph
                    .node_weights()
                    .map(|a| if a.protons == 0 { 85 } else { a.protons })
                    .collect::<Vec<_>>();
                let edges = graph
                    .edge_references()
                    .map(|e| {
                        [
                            e.source().index() as u16,
                            e.target().index() as u16,
                            e.weight().bond_count().floor() as u16,
                        ]
                    })
                    .collect::<Vec<_>>();
                let locs = coordgen::gen_coords(&atoms, &edges).unwrap();
                let mut out = String::new();
                let min_x = locs.iter().map(|p| p.0).min_by(f32::total_cmp).unwrap_or(0.0);
                let min_y = locs.iter().map(|p| p.1).min_by(f32::total_cmp).unwrap_or(0.0);
                let max_x = locs.iter().map(|p| p.0).max_by(f32::total_cmp).unwrap_or(0.0);
                let max_y = locs.iter().map(|p| p.1).max_by(f32::total_cmp).unwrap_or(0.0);
                let diff_x = max_x - min_x + 40.0;
                let diff_y = max_y - min_y + 40.0;
                let max_axis =
                    std::cmp::max_by(diff_x, diff_y, f32::total_cmp);
                let mut add_x = -min_x + 20.0;
                let mut add_y = -min_y + 20.0;
                if diff_y > diff_x {
                    add_x += (diff_y - diff_x) / 2.0;
                } else {
                    add_y += (diff_x - diff_y) / 2.0;
                }
                let locs = locs
                    .iter()
                    .map(|(x, y)| {
                        (
                            (x + add_x).floor() as u16,
                            (y + add_y).floor() as u16,
                        )
                    })
                    .collect::<Vec<_>>();

                for edge in graph.edge_references() {
                    let (x1, y1) = locs[edge.source().index()];
                    let (x2, y2) = locs[edge.target().index()];
                    match edge.weight() {
                        Bond::Non => {},
                        Bond::Single | Bond::Left | Bond::Right => out += &format!("  <line x1=\"{x1}\" y1=\"{y1}\" x2=\"{x2}\" y2=\"{y2}\" style=\"stroke:black;stroke-width:2\"/>\n"),
                        Bond::Double => out += &format!("  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:black;stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:black;stroke-width:2\" />\n", x1 - 2, y1 - 2, x2 - 2, y2 - 2, x1 + 2, y1 + 2, x2 + 2, y2 + 2),
                        Bond::Triple => out += &format!("  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:black;stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:black;stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:black;stroke-width:2\" />\n", x1 - 4, y1 - 4, x2 - 4, y2 - 4, x1 + 4, y1 + 4, x2 + 4, y2 + 4, x1, y1, x2, y2),
                        Bond::Quad => out += &format!("  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:black;stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:black;stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:black;stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:black;stroke-width:2\" />\n", x1 - 6, y1 - 6, x2 - 6, y2 - 6, x1 - 2, y1 - 2, x2 - 2, y2 - 2, x1 + 2, y1 + 2, x2 + 2, y2 + 2, x1 + 6, y1 + 6, x2 + 6, y2 + 6),
                        Bond::Aromatic => out += &format!("  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:black;stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke-dasharray=\"10,10\" style=\"stroke:black;stroke-width:2\" />\n", x1 - 2, y1 - 2, x2 - 2, y2 - 2, x1 + 2, y1 + 2, x2 + 2, y2 + 2),
                    }
                }

                if !edges.is_empty() {
                    out.push('\n')
                }

                for (atom, (cx, cy)) in graph.node_weights().zip(&locs) {
                    let r = match atom.protons {
                        0 => 30,
                        1 | 2 => 8,
                        3..=10 => 14,
                        11..=18 => 18,
                        19..=36 => 22,
                        37..=54 => 26,
                        55.. => 39,
                    };
                    out += &format!(
                        "  <circle r=\"{r}\" cx=\"{cx}\" cy=\"{cy}\" fill=\"{}\" />\n  <text x=\"{}\" y=\"{}\" font-size=\"8\" fill=\"#333\">{atom}</text>\n",
                        svg_atom_color(atom.protons),
                        cx - 3,
                        cy + 3,
                    );
                }
                println!("<svg width=\"{max_axis}\" height=\"{max_axis}\">\n{out}</svg>");
            }
            OutputType::Png => todo!(),
        },
        Err(err) => eprintln!("{err}"),
    }
}
