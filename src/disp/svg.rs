use crate::atom_info::*;
use crate::core::*;
use crate::prelude::DataValueMap;
use fmtastic::*;
use petgraph::visit::*;
use petgraph::Undirected;
use std::fmt::{self, Debug, Display};

#[cfg(feature = "resvg")]
use resvg::*;

pub const SVG_SUPPRESSED_R: &str = "#407F00";
pub const SVG_SUPPRESSED_H: &str = "#577478";
pub const SVG_BOND_COLOR: &str = "#777";

pub fn svg_atom_color(num: u8) -> &'static str {
    match num {
        0 => "#7FFF00",
        1 => "#ADD8E6",
        3 => "#9400D3",
        4 => "#8B008B",
        5 => "#CC7000",
        6 => "#606060",
        7 => "#00BFFF",
        8 => "#CC0000",
        9 => "#00C000",
        11 => "#7900AD",
        12 => "#750075",
        14 => "#663800",
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

fn atom_radius(protons: u8) -> u8 {
    match protons {
        0 => 25,
        1 | 2 => 8,
        3..=10 => 11,
        11..=18 => 14,
        19..=36 => 16,
        37..=54 => 19,
        55.. => 23,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FormatMode {
    /// Normal mode, gives a diagram
    Normal,
    /// Legacy mode, draws atoms as circles
    Legacy,
    /// Legacy mode, with hydrogens explicitly drawn
    LegacyH,
}
impl FormatMode {
    pub const fn is_legacy(&self) -> bool {
        matches!(self, Self::Legacy | Self::LegacyH)
    }
}

pub fn fmt_as_svg<G>(graph: G) -> SvgFormatter<G>
where
    G: Data<NodeWeight = Atom, EdgeWeight = Bond>
        + GraphProp<EdgeType = Undirected>
        + IntoNodeReferences
        + IntoEdgeReferences
        + NodeCompactIndexable,
{
    SvgFormatter {
        graph,
        mode: FormatMode::Normal,
    }
}

#[cfg(feature = "resvg")]
lazy_static::lazy_static! {
    static ref OPTS: usvg::Options = usvg::Options {font_family: "DejaVu Sans Mono".to_string(), ..Default::default()};
    static ref FONTS: usvg::fontdb::Database = {
        use usvg::fontdb::*;
        let mut db = Database::new();
        db.load_font_source(Source::Binary(std::sync::Arc::new(include_bytes!("../../data/DejaVuSansMono.ttf"))));
        db
    };
}

#[derive(Debug, Clone, Copy)]
pub struct SvgFormatter<G> {
    pub graph: G,
    pub mode: FormatMode,
}
impl<G> SvgFormatter<G> {
    pub fn new(graph: G) -> Self {
        Self {
            graph,
            mode: FormatMode::Normal,
        }
    }
}
#[cfg(feature = "resvg")]
impl<G> SvgFormatter<G>
where
    G: DataValueMap<NodeWeight = Atom, EdgeWeight = Bond>
        + GraphProp<EdgeType = Undirected>
        + IntoNodeReferences
        + IntoEdgeReferences
        + IntoEdges
        + NodeCompactIndexable,
{
    pub fn render_to_bytes(&self, bytes: &mut [u8], w: u32, h: u32) {
        let s = self.to_string();
        let tree = usvg::Tree::from_str(&s, &OPTS, &FONTS).unwrap();
        let mut pm = tiny_skia::PixmapMut::from_bytes(bytes, w, h).unwrap();
        resvg::render(&tree, Default::default(), &mut pm);
    }
    pub fn render(&self, size: Option<(u32, u32)>) -> tiny_skia::Pixmap {
        let s = self.to_string();
        let tree = usvg::Tree::from_str(&s, &OPTS, &FONTS).unwrap();
        let (w, h) = size.unwrap_or_else(|| tree.size().to_int_size().dimensions());
        let mut pm = tiny_skia::Pixmap::new(w, h).unwrap();
        resvg::render(&tree, Default::default(), &mut pm.as_mut());
        pm
    }
}
#[cfg(feature = "coordgen")]
impl<G> Display for SvgFormatter<G>
where
    G: DataValueMap<NodeWeight = Atom, EdgeWeight = Bond>
        + GraphProp<EdgeType = Undirected>
        + IntoNodeReferences
        + IntoEdgeReferences
        + IntoEdges
        + NodeCompactIndexable,
{
    #[allow(clippy::write_with_newline)]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut atoms = self
            .graph
            .node_references()
            .map(|a| {
                if a.weight().protons == 0 {
                    85
                } else {
                    a.weight().protons
                }
            })
            .collect::<Vec<_>>();
        let mut edges = self
            .graph
            .edge_references()
            .map(|e| {
                [
                    self.graph.to_index(e.source()) as u16,
                    self.graph.to_index(e.target()) as u16,
                    e.weight().bond_count().floor() as u16,
                ]
            })
            .collect::<Vec<_>>();
        for atom in self.graph.node_references() {
            let idx = self.graph.to_index(atom.id());
            let data = atom.weight().data;
            if self.mode == FormatMode::LegacyH {
                for _ in 0..data.hydrogen() {
                    edges.push([idx as _, atoms.len() as _, 1]);
                    atoms.push(1);
                }
            }
            for _ in 0..data.unknown() {
                edges.push([idx as _, atoms.len() as _, 1]);
                atoms.push(85);
            }
        }
        let locs = coordgen::gen_coords(&atoms, &edges).unwrap();
        let mut out = String::new();
        let (min_x, min_y, max_x, max_y) = if atoms.is_empty() {
            (0.0, 0.0, 0.0, 0.0)
        } else {
            locs.iter().zip(&atoms).fold(
                (
                    f32::INFINITY,
                    f32::INFINITY,
                    f32::NEG_INFINITY,
                    f32::NEG_INFINITY,
                ),
                |(ix, iy, ax, ay), (l, a)| {
                    let r = if self.mode.is_legacy() {
                        atom_radius(*a) as f32
                    } else {
                        0.0
                    };
                    (
                        ix.min(l.0 - r),
                        iy.min(l.1 - r),
                        ax.max(l.0 + r),
                        ay.max(l.1 + r),
                    )
                },
            )
        };
        let diff_x = max_x - min_x + 40.0;
        let diff_y = max_y - min_y + 40.0;
        let max_axis = std::cmp::max_by(diff_x, diff_y, f32::total_cmp);
        writeln!(f, "<svg width=\"{max_axis}px\" height=\"{max_axis}px\" xmlns=\"http://www.w3.org/2000/svg\">")?;
        let mut add_x = -min_x + 20.0;
        let mut add_y = -min_y + 20.0;
        if diff_y > diff_x {
            add_x += (diff_y - diff_x) / 2.0;
        } else {
            add_y += (diff_x - diff_y) / 2.0;
        }

        for edge in self.graph.edge_references() {
            let [ix1, ix2] = std::cmp::minmax(
                self.graph.to_index(edge.source()),
                self.graph.to_index(edge.target()),
            );
            let (ox1, oy1) = locs[ix1];
            let (ox2, oy2) = locs[ix2];
            let (mut x1, mut y1) = (ox1, oy1);
            let (mut x2, mut y2) = (ox2, oy2);
            let (mut dx, mut dy) = (x2 - x1, y2 - y1);
            let mag = (dx * dx + dy * dy).sqrt();
            x1 += add_x;
            y1 += add_y;
            x2 += add_x;
            y2 += add_y;
            dx /= mag;
            dy /= mag;
            if self.mode == FormatMode::Normal {
                if atoms[ix1] != 6 || self.graph.node_weight(edge.source()).unwrap().isotope != 0 {
                    x1 += dx * 12.0;
                    y1 += dy * 12.0;
                }
                if atoms[ix2] != 6 || self.graph.node_weight(edge.target()).unwrap().isotope != 0 {
                    x2 -= dx * 12.0;
                    y2 -= dy * 12.0;
                }
            }
            let (mut rdx, mut rdy) = (-dy * 2.75, dx * 2.75);
            match *edge.weight() {
                Bond::Non => {},
                Bond::Single | Bond::Left | Bond::Right => write!(f, "  <line x1=\"{x1}\" y1=\"{y1}\" x2=\"{x2}\" y2=\"{y2}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\"/>\n")?,
                Bond::Double | Bond::DoubleE | Bond::DoubleZ => write!(f,
                    "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n",
                    x1 - rdx, y1 - rdy,
                    x2 - rdx, y2 - rdy,
                    x1 + rdx, y1 + rdy,
                    x2 + rdx, y2 + rdy)?,
                Bond::Triple => write!(f,
                    "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n",
                    x1 - 2.0 * rdx, y1 - 2.0 * rdy,
                    x2 - 2.0 * rdx, y2 - 2.0 * rdy,
                    x1 + 2.0 * rdx, y1 + 2.0 * rdy,
                    x2 + 2.0 * rdx, y2 + 2.0 * rdy,
                    x1, y1, x2, y2)?,
                Bond::Quad => write!(f,
                    "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n",
                    x1 - 3.0 * rdx, y1 - 3.0 * rdy,
                    x2 - 3.0 * rdx, y2 - 3.0 * rdy,
                    x1 - rdx, y1 - rdy,
                    x2 - rdx, y2 - rdy,
                    x1 + rdx, y1 + rdy,
                    x2 + rdx, y2 + rdy,
                    x1 + 3.0 * rdx, y1 + 3.0 * rdy,
                    x2 + 3.0 * rdx, y2 + 3.0 * rdy)?,
                Bond::Aromatic => {
                    let mut flip = 0.0;
                    let id1 = edge.source();
                    let id2 = edge.target();
                    for e in self.graph.edges(id1).chain(self.graph.edges(id2)) {
                        if *e.weight() != Bond::Aromatic {
                            continue;
                        }
                        if [id1, id2].contains(&e.target()) {
                            continue;
                        }
                        let (x1, y1) = locs[self.graph.to_index(e.source())];
                        let (x2, y2) = locs[self.graph.to_index(e.target())];
                        let ndx = x2 - x1;
                        let ndy = y2 - y1;
                        let dot = dx * ndy - dy * ndx;
                        flip += dot;
                    }
                    if flip > -0.0001 {
                        rdx = -rdx;
                        rdy = -rdy;
                    }
                    write!(f,
                        "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\"  stroke-dasharray=\"10,10\"/>\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n",
                        x1 - rdx, y1 - rdy,
                        x2 - rdx, y2 - rdy,
                        x1 + rdx, y1 + rdy,
                        x2 + rdx, y2 + rdy
                    )?
                },
                _ => panic!("invalid bond!")
            }
        }

        if !edges.is_empty() {
            out.push('\n')
        }

        let mut idx = self.graph.node_count();

        for (atom, &(mut x1, mut y1)) in self.graph.node_references().zip(&locs) {
            let atom = atom.weight();
            let data = atom.data;
            x1 += add_x;
            y1 += add_y;
            if self.mode == FormatMode::LegacyH {
                for i in 0..data.hydrogen() {
                    let (mut x2, mut y2) = locs[idx + (i as usize)];
                    x2 += add_x;
                    y2 += add_y;
                    let r = atom_radius(1);
                    write!(f, "  <line x1=\"{x1}\" y1=\"{y1}\" x2=\"{x2}\" y2=\"{y2}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\"/>\n  <circle r=\"{r}\" cx=\"{x2}\" cy=\"{y2}\" fill=\"{}\" />\n  <text x=\"{x2}\" y=\"{y2}\" font-size=\"{r}\" text-anchor=\"middle\" alignment-baseline=\"middle\" fill=\"#444\">H</text>\n", SVG_SUPPRESSED_H)?;
                }
                idx += data.hydrogen() as usize;
            }
            for i in 0..data.unknown() {
                let (mut x2, mut y2) = locs[idx + (i as usize)];
                x2 += add_x;
                y2 += add_y;
                if self.mode.is_legacy() {
                    let r = atom_radius(0);
                    write!(f, "  <line x1=\"{x1}\" y1=\"{y1}\" x2=\"{x2}\" y2=\"{y2}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\"/>\n  <circle r=\"{r}\" cx=\"{x2}\" cy=\"{y2}\" fill=\"{}\" />\n  <text x=\"{x2}\" y=\"{y2}\" font-size=\"{r}\" text-anchor=\"middle\" alignment-baseline=\"middle\" fill=\"#444\">R</text>\n", SVG_SUPPRESSED_R)?;
                } else {
                    let dx = x2 - x1;
                    let dy = y2 - y1;
                    let mag = (dx * dx + dy * dy).sqrt();
                    let x3 = x2 - dx / mag * 12.0;
                    let y3 = y2 - dy / mag * 12.0;
                    write!(f, "  <line x1=\"{x1}\" y1=\"{y1}\" x2=\"{x3}\" y2=\"{y3}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\"/>\n  <text x=\"{x2}\" y=\"{y2}\" font-size=\"20\" text-anchor=\"middle\" alignment-baseline=\"middle\" fill=\"{SVG_SUPPRESSED_R}\">R</text>\n")?;
                }
            }
            idx += data.unknown() as usize;
        }

        if idx != self.graph.node_count() {
            out.push('\n');
        }

        for (aref, &(mut cx, mut cy)) in self.graph.node_references().zip(&locs) {
            let atom = aref.weight();
            cx += add_x;
            cy += add_y;
            if self.mode.is_legacy() {
                let r = atom_radius(atom.protons);
                write!(f, "  <circle r=\"{r}\" cx=\"{cx}\" cy=\"{cy}\" fill=\"{}\" />\n  <text x=\"{cx}\" y=\"{cy}\" font-size=\"{r}\" text-anchor=\"middle\" alignment-baseline=\"middle\" fill=\"#444\">{atom}</text>\n", svg_atom_color(atom.protons))?;
            } else {
                let left = [8, 9, 16, 17, 34, 35, 52, 53, 84, 85].contains(&atom.protons);
                let color = svg_atom_color(atom.protons);
                let mut tx = cx;
                let w = ATOM_DATA[atom.protons as usize].sym.len()
                    + match atom.data.hydrogen() {
                        0 => 0,
                        n => n.ilog10() as usize + 1,
                    }
                    + match atom.isotope {
                        0 => 0,
                        n => n.ilog10() as usize + 1,
                    };
                tx -= w as f32 * 5.0;
                if let Some((n, (mut dx, mut dy))) = self
                    .graph
                    .neighbors(aref.id())
                    .map(|i| (1usize, locs[self.graph.to_index(i)]))
                    .reduce(|(a, (ax, ay)), (b, (bx, by))| (a + b, (ax + bx, ay + by)))
                {
                    dx /= n as f32;
                    dy /= n as f32;
                    dx -= cx - add_x;
                    dy -= cy - add_y;
                    let mag = (dx * dx + dy * dy).sqrt();
                    dx /= mag;
                    dy /= mag;
                    if atom.protons != 6 || atom.isotope != 0 {
                        write!(f, "  <text x=\"{tx}\" y=\"{cy}\" font-size=\"15\" alignment-baseline=\"middle\" fill=\"{color}\">")?;
                        if dx > 0.0 && left {
                            match atom.data.hydrogen() {
                                0 => {}
                                1 => f.write_str("H")?,
                                n => write!(f, "H{}", Subscript(n))?,
                            }
                        }
                        if atom.protons == 0 {
                            f.write_str(match atom.isotope {
                                0xFFFE => "A",
                                0xFFFC => "Q",
                                0x0100 => "X",
                                0x042B => "M",
                                _ => "*",
                            })?;
                        } else {
                            if atom.isotope != 0 {
                                write!(f, "{}", Superscript(atom.isotope))?;
                            }
                            f.write_str(ATOM_DATA[atom.protons as usize].sym)?;
                        }
                        if dx < 0.0 || !left {
                            match atom.data.hydrogen() {
                                0 => {}
                                1 => f.write_str("H")?,
                                n => write!(f, "H{}", Subscript(n))?,
                            }
                        }
                        if atom.protons != 6 {
                            match atom.charge {
                                0 => {}
                                1 => f.write_str("⁺")?,
                                -1 => f.write_str("⁻")?,
                                _ => write!(f, "{:+}", Superscript(atom.charge))?,
                            }
                        }
                        f.write_str("</text>\n")?;
                    }
                    if atom.protons == 6 && atom.charge != 0 {
                        write!(f, "  <text x=\"{}\" y=\"{}\" font-size=\"15\" alignment-baseline=\"middle\" fill=\"{color}\">", cx - dx, cy - dy)?;
                        match atom.charge {
                            0 => {}
                            1 => f.write_str("⁺")?,
                            -1 => f.write_str("⁻")?,
                            _ => write!(f, "{:+}", Superscript(atom.charge))?,
                        }
                        f.write_str("</text>\n")?;
                    }
                } else {
                    write!(f, "  <text x=\"{tx}\" y=\"{cy}\" font-size=\"15\" alignment-baseline=\"middle\" fill=\"{color}\">")?;
                    if left {
                        match atom.data.hydrogen() {
                            0 => {}
                            1 => f.write_str("H")?,
                            n => write!(f, "H{}", Subscript(n))?,
                        }
                    }
                    if atom.protons == 0 {
                        f.write_str(match atom.isotope {
                            0xFFFE => "A",
                            0xFFFC => "Q",
                            0x0100 => "X",
                            0x042B => "M",
                            _ => "*",
                        })?;
                    } else {
                        if atom.isotope != 0 {
                            write!(f, "{}", Superscript(atom.isotope))?;
                        }
                        f.write_str(ATOM_DATA[atom.protons as usize].sym)?;
                    }
                    if !left {
                        match atom.data.hydrogen() {
                            0 => {}
                            1 => f.write_str("H")?,
                            n => write!(f, "H{}", Subscript(n))?,
                        }
                    }
                    match atom.charge {
                        0 => {}
                        1 => f.write_str("⁺")?,
                        -1 => f.write_str("⁻")?,
                        _ => write!(f, "{:+}", Superscript(atom.charge))?,
                    }
                    f.write_str("</text>\n")?;
                }
            }
        }

        f.write_str("</svg>")
    }
}
