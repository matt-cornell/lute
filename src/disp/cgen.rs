pub use super::color::*;
use crate::core::*;
use petgraph::visit::*;
use petgraph::Undirected;
use std::fmt::{self, Display, Formatter};

#[cfg(feature = "resvg")]
use resvg::*;

pub fn fmt_with_cg<G>(graph: G) -> CgFormatter<G>
where
    G: Data<NodeWeight = Atom, EdgeWeight = Bond>
        + GraphProp<EdgeType = Undirected>
        + IntoNodeReferences
        + IntoEdgeReferences
        + NodeCompactIndexable,
{
    CgFormatter(graph)
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
pub struct CgFormatter<G>(pub G);
#[cfg(feature = "resvg")]
impl<G> CgFormatter<G>
where
    G: Data<NodeWeight = Atom, EdgeWeight = Bond>
        + IntoNodeReferences
        + IntoEdgeReferences
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
impl<G> Display for CgFormatter<G>
where
    G: Data<NodeWeight = Atom, EdgeWeight = Bond>
        + IntoNodeReferences
        + IntoEdgeReferences
        + NodeCompactIndexable,
{
    #[allow(clippy::write_with_newline)]
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut atoms = self
            .0
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
            .0
            .edge_references()
            .map(|e| {
                [
                    self.0.to_index(e.source()) as u16,
                    self.0.to_index(e.target()) as u16,
                    e.weight().bond_count().floor() as u16,
                ]
            })
            .collect::<Vec<_>>();
        for atom in self.0.node_references() {
            let idx = self.0.to_index(atom.id());
            let data = atom.weight().data;
            for _ in 0..data.hydrogen() {
                edges.push([idx as _, atoms.len() as _, 1]);
                atoms.push(1);
            }
            for _ in 0..data.unknown() {
                edges.push([idx as _, atoms.len() as _, 1]);
                atoms.push(85);
            }
        }
        let locs = coordgen::gen_coords(&atoms, &edges).unwrap();
        let mut out = String::new();
        let min_x = locs
            .iter()
            .zip(&atoms)
            .map(|(p, &a)| p.0 - atom_radius(a) as f32)
            .min_by(f32::total_cmp)
            .unwrap_or(0.0);
        let min_y = locs
            .iter()
            .zip(&atoms)
            .map(|(p, &a)| p.1 - atom_radius(a) as f32)
            .min_by(f32::total_cmp)
            .unwrap_or(0.0);
        let max_x = locs
            .iter()
            .zip(&atoms)
            .map(|(p, &a)| p.0 + atom_radius(a) as f32)
            .max_by(f32::total_cmp)
            .unwrap_or(0.0);
        let max_y = locs
            .iter()
            .zip(&atoms)
            .map(|(p, &a)| p.1 + atom_radius(a) as f32)
            .max_by(f32::total_cmp)
            .unwrap_or(0.0);
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
        let locs = locs
            .iter()
            .map(|(x, y)| ((x + add_x).floor() as i16, (y + add_y).floor() as i16))
            .collect::<Vec<_>>();

        for edge in self.0.edge_references() {
            let (x1, y1) = locs[self.0.to_index(edge.source())];
            let (x2, y2) = locs[self.0.to_index(edge.target())];
            let (dx, dy) = (x2 - x1, y2 - y1);
            let (dx, dy) = (-dy as f64, dx as f64);
            let mag = (dx * dx + dy * dy).sqrt() / 3.0;
            let (dx, dy) = ((dx / mag) as i16, (dy / mag) as i16);
            match *edge.weight() {
                Bond::Non => {},
                Bond::Single | Bond::Left | Bond::Right => write!(f, "  <line x1=\"{x1}\" y1=\"{y1}\" x2=\"{x2}\" y2=\"{y2}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\"/>\n")?,
                Bond::Double | Bond::DoubleE | Bond::DoubleZ => write!(f,
                    "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n",
                    x1 - dx, y1 - dy,
                    x2 - dx, y2 - dy,
                    x1 + dx, y1 + dy,
                    x2 + dx, y2 + dy)?,
                Bond::Triple => write!(f,
                    "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n",
                    x1 - 2 * dx, y1 - 2 * dy,
                    x2 - 2 * dx, y2 - 2 * dy,
                    x1 + 2 * dx, y1 + 2 * dy,
                    x2 + 2 * dx, y2 + 2 * dy,
                    x1, y1, x2, y2)?,
                Bond::Quad => write!(f,
                    "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n",
                    x1 - 3 * dx, y1 - 3 * dy,
                    x2 - 3 * dx, y2 - 3 * dy,
                    x1 - dx, y1 - dy,
                    x2 - dx, y2 - dy,
                    x1 + dx, y1 + dy,
                    x2 + dx, y2 + dy,
                    x1 + 3 * dx, y1 + 3 * dy,
                    x2 + 3 * dx, y2 + 3 * dy)?,
                Bond::Aromatic => write!(f,
                    "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\"  stroke-dasharray=\"10,10\"/>\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n",
                    x1 - dx, y1 - dy,
                    x2 - dx, y2 - dy,
                    x1 + dx, y1 + dy,
                    x2 + dx, y2 + dy)?,
                _ => panic!("invalid bond!")
            }
        }

        if !edges.is_empty() {
            out.push('\n')
        }

        let mut idx = self.0.node_count();

        for (atom, (x1, y1)) in self.0.node_references().zip(&locs) {
            let atom = atom.weight();
            let data = atom.data;
            for i in 0..data.hydrogen() {
                let (x2, y2) = locs[idx + (i as usize)];
                let r = atom_radius(1);
                write!(f, "  <line x1=\"{x1}\" y1=\"{y1}\" x2=\"{x2}\" y2=\"{y2}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\"/>\n  <circle r=\"{r}\" cx=\"{x2}\" cy=\"{y2}\" fill=\"{}\" />\n  <text x=\"{x2}\" y=\"{y2}\" font-size=\"{r}\" text-anchor=\"middle\" alignment-baseline=\"middle\" fill=\"#444\">H</text>\n", SVG_SUPPRESSED_H)?;
            }
            idx += data.hydrogen() as usize;
            for i in 0..data.unknown() {
                let (x2, y2) = locs[idx + (i as usize)];
                let r = atom_radius(0);
                write!(f, "  <line x1=\"{x1}\" y1=\"{y1}\" x2=\"{x2}\" y2=\"{y2}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\"/>\n  <circle r=\"{r}\" cx=\"{x2}\" cy=\"{y2}\" fill=\"{}\" />\n  <text x=\"{x2}\" y=\"{y2}\" font-size=\"{r}\" text-anchor=\"middle\" alignment-baseline=\"middle\" fill=\"#444\">R</text>\n", SVG_SUPPRESSED_R)?;
            }
            idx += data.unknown() as usize;
        }

        if idx != self.0.node_count() {
            out.push('\n');
        }

        for (atom, (cx, cy)) in self.0.node_references().zip(&locs) {
            let atom = atom.weight();
            let r = atom_radius(atom.protons);
            write!(f, "  <circle r=\"{r}\" cx=\"{cx}\" cy=\"{cy}\" fill=\"{}\" />\n  <text x=\"{cx}\" y=\"{cy}\" font-size=\"{r}\" text-anchor=\"middle\" alignment-baseline=\"middle\" fill=\"#444\">{atom}</text>\n", svg_atom_color(atom.protons))?;
        }

        f.write_str("</svg>")
    }
}
