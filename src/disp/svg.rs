use crate::atom_info::*;
use crate::core::*;
use coordgen_rs::prelude::*;
use coordgen_rs::sketcher::PointF;
use petgraph::visit::*;
use petgraph::Undirected;
use std::fmt::{self, Display, Formatter};


#[cfg(feature = "resvg")]
use resvg::*;

pub const SVG_SUPPRESSED_R: &str = "#407F00";
pub const SVG_SUPPRESSED_H: &str = "#577478";
pub const SVG_BOND_COLOR: &str = "#444";
pub fn svg_atom_color(num: u8) -> &'static str {
    match num {
        0 => "#7FFF00",
        1 => "#ADD8E6",
        3 => "#9400D3",
        4 => "#8B008B",
        5 => "#CC7000",
        6 => "#0F0F0F",
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

pub fn fmt_as_svg<G>(graph: G) -> SvgFormatter<G>
where
    G: Data<NodeWeight = Atom, EdgeWeight = Bond>
        + GraphProp<EdgeType = Undirected>
        + IntoNodeReferences
        + IntoEdgeReferences
        + NodeCompactIndexable,
{
    SvgFormatter(graph)
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
pub struct SvgFormatter<G>(pub G);
#[cfg(feature = "resvg")]
impl<G> SvgFormatter<G>
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
impl<G> Display for SvgFormatter<G>
where
    G: Data<NodeWeight = Atom, EdgeWeight = Bond>
        + IntoNodeReferences
        + IntoEdgeReferences
        + NodeCompactIndexable,
{
    #[allow(clippy::write_with_newline)]
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let intern = Unsync::new();
        let mut builder = Builder::new(&intern);
        let mut atoms = vec![None; self.0.node_count()];
        for a in self.0.node_references() {
            let p = if a.weight().protons == 0 {
                85
            } else {
                a.weight().protons
            };
            let r = builder.add_atom(p);
            r.borrow_mut().charge = a.weight().charge;
            atoms[self.0.to_index(a.id())] = Some(r);
        }
        for e in self.0.edge_references() {
            builder.add_bond(
                atoms[self.0.to_index(e.source())].unwrap(),
                atoms[self.0.to_index(e.source())].unwrap(),
                e.weight().bond_count().floor() as _,
            );
        }
        for atom in self.0.node_references() {
            let idx = self.0.to_index(atom.id());
            let data = atom.weight().data;
            for _ in 0..data.unknown() {
                let h = builder.add_atom(1);
                builder.add_bond(atoms[idx as usize].unwrap(), h, 85);
            }
        }
        let mol = builder.finish();
        let mut sketcher = Sketcher::new(&intern);
        sketcher.generate(mol);
        let (min_x, min_y, max_x, max_y) = sketcher.atoms().iter().map(|a| {
            let PointF(x, y) = a.borrow().coordinates;
            (x - 20.0, x + 20.0, y - 20.0, y + 20.0)
        }).reduce(|(ix1, iy1, ax1, ay1), (ix2, iy2, ax2, ay2)| {
            (ix1.min(ix2), iy1.min(iy2), ax1.max(ax2), ay1.max(ay2))
        }).unwrap_or((0.0, 0.0, 0.0, 0.0));
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

        for edge in self.0.edge_references() {
            let a1 = atoms[self.0.to_index(edge.source())].unwrap().borrow();
            let a2 = atoms[self.0.to_index(edge.target())].unwrap().borrow();
            let mut c1 = a1.coordinates;
            let mut c2 = a2.coordinates;
            let mut d = c2 - c1;
            d.normalize();
            if a1.atom_number != 6 {
                c1 += d * 10.0;
            }
            if a2.atom_number != 6 {
                c2 -= d * 10.0;
            }
            let PointF(mut x1, mut y1) = c1;
            let PointF(mut x2, mut y2) = c2;
            let PointF(mut dy, dx) = d * 0.3;
            dy = -dy;
            x1 += add_x;
            y1 += add_y;
            x2 += add_x;
            y2 += add_y;
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
                    x1 - 2.0 * dx, y1 - 2.0 * dy,
                    x2 - 2.0 * dx, y2 - 2.0 * dy,
                    x1 + 2.0 * dx, y1 + 2.0 * dy,
                    x2 + 2.0 * dx, y2 + 2.0 * dy,
                    x1, y1, x2, y2)?,
                Bond::Quad => write!(f,
                    "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n",
                    x1 - 3.0 * dx, y1 - 3.0 * dy,
                    x2 - 3.0 * dx, y2 - 3.0 * dy,
                    x1 - dx, y1 - dy,
                    x2 - dx, y2 - dy,
                    x1 + dx, y1 + dy,
                    x2 + dx, y2 + dy,
                    x1 + 3.0 * dx, y1 + 3.0 * dy,
                    x2 + 3.0 * dx, y2 + 3.0 * dy)?,
                Bond::Aromatic => write!(f,
                    "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\"  stroke-dasharray=\"10,10\"/>\n  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\" />\n",
                    x1 - dx, y1 - dy,
                    x2 - dx, y2 - dy,
                    x1 + dx, y1 + dy,
                    x2 + dx, y2 + dy)?,
                _ => panic!("invalid bond!")
            }
        }

        let mut idx = self.0.node_count();

        for a in self.0.node_references() {
            let atom = a.weight();
            let data = atom.data;
            let PointF(x1, y1) = atoms[self.0.to_index(a.id())].unwrap().borrow().coordinates;
            write!(f, "  <text x=\"{x1}\" y=\"{y1}\" font-size=\"10\" text-anchor=\"middle\" alignment-baseline=\"middle\" fill=\"{}\">{atom}</text>\n", svg_atom_color(atom.protons))?;
            for i in 0..data.unknown() {
                let PointF(x2, y2) = atoms[idx + (i as usize)].unwrap().borrow().coordinates;
                write!(f, "  <line x1=\"{x1}\" y1=\"{y1}\" x2=\"{x2}\" y2=\"{y2}\" style=\"stroke:{SVG_BOND_COLOR};stroke-width:2\"/>\n  <text x=\"{x2}\" y=\"{y2}\" font-size=\"10\" text-anchor=\"middle\" alignment-baseline=\"middle\" fill=\"{SVG_SUPPRESSED_R}\">R</text>\n")?;
            }
            idx += data.unknown() as usize;
        }

        f.write_str("</svg>")
    }
}
