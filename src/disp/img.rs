use crate::atom_info::*;
use crate::molecule::*;
use fimg::Image;
use petgraph::visit::*;
use petgraph::Undirected;
use std::num::NonZeroU32;

const RAW_FONT_DATA: &[u8] = include_bytes!("../../data/LiberationMono-Regular.ttf");

lazy_static::lazy_static! {
    pub static ref FONT_DATA: fontdue::Font = fontdue::Font::from_bytes(RAW_FONT_DATA, Default::default()).expect("font loading should not fail!");
}

pub const IMG_SUPPRESSED_R: u32 = 0x407F00;
pub const IMG_SUPPRESSED_H: u32 = 0x577478;
pub const IMG_BOND_COLOR: u32 = 0x404040;
pub const IMG_ATOM_TEXT: u32 = 0x4F4F4F;

pub fn img_atom_color(num: u8) -> u32 {
    match num {
        0 => 0x7FFF00,
        1 => 0xADD8E6,
        3 => 0x9400D3,
        4 => 0x8B008B,
        5 => 0xCC7000,
        6 => 0x0F0F0F,
        7 => 0x00BFFF,
        8 => 0xCC0000,
        9 => 0x00C000,
        11 => 0x7900AD,
        12 => 0x750075,
        14 => 0x663800,
        15 => 0x8B0000,
        16 => 0xDAA520,
        17 => 0x00A000,
        19 => 0x60008A,
        20 => 0x610061,
        35 => 0x008000,

        _ => match ATOM_DATA[num as usize].group {
            ElemGroup::Nonmet => 0x20B2AA,
            ElemGroup::Noble => 0xFF69B4,
            ElemGroup::Trans => 0x778899,
            ElemGroup::Poor => 0x2F4F4F,
            ElemGroup::RarEar => 0xFF1943,
            ElemGroup::Metoid => 0x556B2F,
            ElemGroup::Alkali => 0x470066,
            ElemGroup::AlkEar => 0x610061,
            ElemGroup::Halogn => 0x006000,
        },
    }
}

pub fn atom_radius(protons: u8) -> u8 {
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

/// Convenience wrapper around `make_img` to allocate a `Vec<u8>`.
pub fn make_img_vec<G>(graph: G) -> Image<Vec<u8>, 4>
where
    G: Data<NodeWeight = Atom, EdgeWeight = Bond>
        + GraphProp<EdgeType = Undirected>
        + GraphRef
        + IntoNodeReferences
        + IntoEdgeReferences
        + NodeCompactIndexable,
{
    make_img(graph, |size| vec![0u8; (size * size * 4) as usize])
}

/// Create an image. Take a callback to generate the buffer for efficiency. `gen` takes the length
/// of the square screen and returns a buffer with a length equal to `4 * len * len`.
pub fn make_img<G, B, F>(graph: G, gen: F) -> Image<B, 4>
where
    B: AsRef<[u8]> + AsMut<[u8]>,
    F: FnOnce(u32) -> B,
    G: Data<NodeWeight = Atom, EdgeWeight = Bond>
        + GraphProp<EdgeType = Undirected>
        + GraphRef
        + IntoNodeReferences
        + IntoEdgeReferences
        + NodeCompactIndexable,
{
    let mut atoms = graph
        .node_references()
        .map(|a| {
            if a.weight().protons == 0 {
                85
            } else {
                a.weight().protons
            }
        })
        .collect::<Vec<_>>();
    let mut edges = graph
        .edge_references()
        .map(|e| {
            [
                graph.to_index(e.source()) as u16,
                graph.to_index(e.target()) as u16,
                e.weight().bond_count().floor() as u16,
            ]
        })
        .collect::<Vec<_>>();
    for atom in graph.node_references() {
        let idx = graph.to_index(atom.id());
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
    let mut add_x = -min_x + 20.0;
    let mut add_y = -min_y + 20.0;
    if diff_y > diff_x {
        add_x += (diff_y - diff_x) / 2.0;
    } else {
        add_y += (diff_x - diff_y) / 2.0;
    }
    let locs = locs
        .iter()
        .map(|(x, y)| ((x + add_x).floor() as i32, (y + add_y).floor() as i32))
        .collect::<Vec<_>>();

    let ma = max_axis.ceil() as u32;
    let mut img = unsafe {
        assert_ne!(ma, 0);
        let buf = gen(ma);
        assert_eq!(buf.as_ref().len(), (ma * ma * 4) as usize);
        Image::new(NonZeroU32::new_unchecked(ma), NonZeroU32::new_unchecked(ma), buf)
    };

    for edge in graph.edge_references() {
        let (x1, y1) = locs[graph.to_index(edge.source())];
        let (x2, y2) = locs[graph.to_index(edge.target())];
        let (dx, dy) = (x2 - x1, y2 - y1);
        let (dx, dy) = (-dy as f64, dx as f64);
        let mag = (dx * dx + dy * dy).sqrt() / 3.0;
        let (dx, dy) = ((dx / mag) as i32, (dy / mag) as i32);
        let mut color = IMG_BOND_COLOR.to_le_bytes();
        match *edge.weight() {
            Bond::Non => {}
            Bond::Single | Bond::Left | Bond::Right => img.line((x1, x2), (y1, y2), color),
            Bond::Double => {
                img.line((x1 - dx, y1 - dy), (x2 - dx, y2 - dy), color);
                img.line((x1 + dx, y1 + dy), (x2 + dx, y2 + dy), color);
            }
            Bond::Triple => {
                img.line((x1, x2), (y1, y2), color);
                img.line(
                    (x1 - 2 * dx, y1 - 2 * dy),
                    (x2 - 2 * dx, y2 - 2 * dy),
                    color,
                );
                img.line(
                    (x1 + 2 * dx, y1 + 2 * dy),
                    (x2 + 2 * dx, y2 + 2 * dy),
                    color,
                );
            }
            Bond::Quad => {
                img.line((x1 - dx, y1 - dy), (x2 - dx, y2 - dy), color);
                img.line((x1 + dx, y1 + dy), (x2 + dx, y2 + dy), color);
                img.line(
                    (x1 - 3 * dx, y1 - 3 * dy),
                    (x2 - 3 * dx, y2 - 3 * dy),
                    color,
                );
                img.line(
                    (x1 + 3 * dx, y1 + 3 * dy),
                    (x2 + 3 * dx, y2 + 3 * dy),
                    color,
                );
            }
            Bond::Aromatic => {
                img.line((x1 - dx, y1 - dy), (x2 - dx, y2 - dy), color);
                color[3] /= 2;
                img.line((x1 + dx, y1 + dy), (x2 + dx, y2 + dy), color);
            }
            _ => panic!("invalid bond!"),
        }
    }

    let mut idx = graph.node_count();

    let tc = IMG_ATOM_TEXT.to_le_bytes();

    for (atom, &(x1, y1)) in graph.node_references().zip(&locs) {
        let atom = atom.weight();
        let data = atom.data;
        let bond = IMG_BOND_COLOR.to_le_bytes();
        let hc = IMG_SUPPRESSED_H.to_le_bytes();
        let rc = IMG_SUPPRESSED_R.to_le_bytes();
        {
            let r = atom_radius(1);
            let metrics = FONT_DATA.metrics('H', r as _);
            for i in 0..data.hydrogen() {
                let (x2, y2) = locs[idx + (i as usize)];
                img.line((x1, y1), (x2, y2), bond);
                img.circle((x2, y2), r as _, hc);
                img.text(
                    (x2 as u32) - (metrics.width as u32) / 2,
                    (y2 as u32) - (metrics.height as u32) / 2,
                    r as _,
                    &FONT_DATA,
                    "H",
                    tc,
                );
            }
        }
        idx += data.hydrogen() as usize;
        {
            let r = atom_radius(0);
            let metrics = FONT_DATA.metrics('R', r as _);
            for i in 0..data.unknown() {
                let (x2, y2) = locs[idx + (i as usize)];
                img.line((x1, y1), (x2, y2), bond);
                img.circle((x2, y2), r as _, rc);
                img.text(
                    (x2 as u32) - (metrics.width as u32) / 2,
                    (y2 as u32) - (metrics.height as u32) / 2,
                    r as _,
                    &FONT_DATA,
                    "R",
                    tc,
                );
            }
        }
        idx += data.unknown() as usize;
    }

    for (atom, &(cx, cy)) in graph.node_references().zip(&locs) {
        let atom = atom.weight();
        let r = atom_radius(atom.protons);

        let metrics = FONT_DATA.metrics('H', r as _);
        let txt = atom.to_string();
        img.text(
            (cx as u32) - ((metrics.width * txt.len()) as u32) / 2,
            (cy as u32) - (metrics.height as u32) / 2,
            r as _,
            &FONT_DATA,
            &txt,
            tc,
        );
    }

    img
}
