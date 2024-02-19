use crate::atom_info::*;
use crate::molecule::*;
use fimg::Image;
use image::*;
use petgraph::visit::*;
use petgraph::Undirected;
use std::num::NonZeroU32;
use std::ops::Deref;

const RAW_FONT_DATA: &[u8] = include_bytes!("../../data/Font.ttf");

lazy_static::lazy_static! {
    pub static ref FONT_DATA: fontdue::Font = fontdue::Font::from_bytes(RAW_FONT_DATA, Default::default()).expect("font loading should not fail!");
}

pub const IMG_SUPPRESSED_R: u32 = 0x407F00FF;
pub const IMG_SUPPRESSED_H: u32 = 0x577478FF;
pub const IMG_BOND_COLOR: u32 = 0x404040FF;
pub const IMG_ATOM_TEXT: u32 = 0x4F4F4FFF;

pub fn img_atom_color(num: u8) -> u32 {
    match num {
        0 => 0x7FFF00FF,
        1 => 0xADD8E6FF,
        3 => 0x9400D3FF,
        4 => 0x8B008BFF,
        5 => 0xCC7000FF,
        6 => 0x0F0F0FFF,
        7 => 0x00BFFFFF,
        8 => 0xCC0000FF,
        9 => 0x00C000FF,
        11 => 0x7900ADFF,
        12 => 0x750075FF,
        14 => 0x663800FF,
        15 => 0x8B0000FF,
        16 => 0xDAA520FF,
        17 => 0x00A000FF,
        19 => 0x60008AFF,
        20 => 0x610061FF,
        35 => 0x008000FF,

        _ => match ATOM_DATA[num as usize].group {
            ElemGroup::Nonmet => 0x20B2AAFF,
            ElemGroup::Noble => 0xFF69B4FF,
            ElemGroup::Trans => 0x778899FF,
            ElemGroup::Poor => 0x2F4F4FF,
            ElemGroup::RarEar => 0xFF1943FF,
            ElemGroup::Metoid => 0x556B2FFF,
            ElemGroup::Alkali => 0x470066FF,
            ElemGroup::AlkEar => 0x610061FF,
            ElemGroup::Halogn => 0x006000FF,
        },
    }
}

pub fn atom_radius(protons: u8) -> u8 {
    match protons {
        0 => 50,
        1 | 2 => 16,
        3..=10 => 22,
        11..=18 => 28,
        19..=36 => 32,
        37..=54 => 38,
        55.. => 36,
    }
}

/// Convenience wrapper around `make_img` to allocate a `Vec<u8>`.
pub fn make_img_vec<G>(graph: G) -> ImageBuffer<Rgba<u8>, Vec<u8>>
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
pub fn make_img<G, B, F>(graph: G, gen: F) -> ImageBuffer<Rgba<u8>, B>
where
    B: AsRef<[u8]> + AsMut<[u8]> + Deref<Target = [u8]>,
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
    let mut locs = coordgen::gen_coords(&atoms, &edges).unwrap();
    let mut min_x = 0.0;
    let mut min_y = 0.0;
    let mut max_x = 0.0;
    let mut max_y = 0.0;
    for ((x, y), a) in locs.iter_mut().zip(&atoms) {
        *x *= 2.0;
        *y *= 2.0;
        let r = atom_radius(*a) as f32;
        if *x - r < min_x {
            min_x = *x - r;
        }
        if *x + r > max_x {
            max_x = *x + r;
        }
        if *y - r < min_y {
            min_y = *y - r;
        }
        if *y + r > max_y {
            max_y = *y + r;
        }
    }

    let diff_x = max_x - min_x + 80.0;
    let diff_y = max_y - min_y + 80.0;
    let max_axis = std::cmp::max_by(diff_x, diff_y, f32::total_cmp);
    let mut add_x = -min_x + 40.0;
    let mut add_y = -min_y + 40.0;
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
        Image::new(
            NonZeroU32::new_unchecked(ma),
            NonZeroU32::new_unchecked(ma),
            buf,
        )
    };

    for edge in graph.edge_references() {
        let (x1, y1) = locs[graph.to_index(edge.source())];
        let (x2, y2) = locs[graph.to_index(edge.target())];
        let (dx, dy) = (x2 - x1, y2 - y1);
        let (dx, dy) = (-dy as f64, dx as f64);
        let mag = (dx * dx + dy * dy).sqrt() / 3.0;
        let (dx, dy) = ((dx / mag) as i32, (dy / mag) as i32);
        let mut color = IMG_BOND_COLOR.to_be_bytes();
        match *edge.weight() {
            Bond::Non => {}
            Bond::Single | Bond::Left | Bond::Right => img.line((x1, y1), (x2, y2), color),
            Bond::Double => {
                img.line((x1 - dx, y1 - dy), (x2 - dx, y2 - dy), color);
                img.line((x1 + dx, y1 + dy), (x2 + dx, y2 + dy), color);
            }
            Bond::Triple => {
                img.line((x1, y1), (x2, y2), color);
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

    let tc = IMG_ATOM_TEXT.to_be_bytes();
    for (atom, &(x1, y1)) in graph.node_references().zip(&locs) {
        let atom = atom.weight();
        let data = atom.data;
        let bond = IMG_BOND_COLOR.to_be_bytes();
        let hc = IMG_SUPPRESSED_H.to_be_bytes();
        let rc = IMG_SUPPRESSED_R.to_be_bytes();
        {
            let r = atom_radius(1);
            let f = r as f32 * 0.9;
            let metrics = FONT_DATA.metrics('H', f);
            for i in 0..data.hydrogen() {
                let (x2, y2) = locs[idx + (i as usize)];
                img.line((x1, y1), (x2, y2), bond);
                img.circle((x2, y2), r as _, hc);
                img.text(
                    (x2 as u32) - (metrics.width as u32) / 2,
                    (y2 as u32) - (metrics.height as u32) * 2 / 3,
                    f,
                    &FONT_DATA,
                    "H",
                    tc,
                );
            }
        }
        idx += data.hydrogen() as usize;
        {
            let r = atom_radius(0);
            let f = r as f32 * 0.9;
            let metrics = FONT_DATA.metrics('R', f);
            for i in 0..data.unknown() {
                let (x2, y2) = locs[idx + (i as usize)];
                img.line((x1, y1), (x2, y2), bond);
                img.circle((x2, y2), r as _, rc);
                img.text(
                    (x2 as u32) - (metrics.width as u32) / 2,
                    (y2 as u32) - (metrics.height as u32) * 2 / 3,
                    f,
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
        let f = r as f32 * 0.9;
        let mut x = 0;
        let mut y = 0;
        let txt = ATOM_DATA[atom.protons as usize].sym;
        for c in txt.chars() {
            let metrics = FONT_DATA.metrics(c, f);
            x += metrics.width;
            if metrics.height > y {
                y = metrics.height;
            }
        }
        img.circle((cx, cy), r as _, img_atom_color(atom.protons).to_be_bytes());
        let top = (cy as u32) - (y as u32) * 2 / 3;
        let lef = (cx as u32) - (x as u32) / 2;
        let rig = (cx as u32) + (x as u32) / 2;

        if atom.isotope != 0 || atom.protons == 0 {
            let f = f / 2.0;
            let txt = atom.isotope.to_string();
            let mut ix = 0;
            let mut iy = 0;
            for c in txt.chars() {
                let metrics = FONT_DATA.metrics(c, f);
                ix += metrics.width;
                if metrics.height > iy {
                    iy = metrics.height;
                }
            }
            img.text(lef - ix as u32 - 2, top, f, &FONT_DATA, &txt, tc);
        }

        if atom.charge != 0 {
            let f = f / 2.0;
            let txt = match atom.charge {
                -1 => "-".to_string(),
                1 => "+".to_string(),
                c => format!("{c:+}"),
            };
            let mut iy = 0;
            for c in txt.chars() {
                let metrics = FONT_DATA.metrics(c, f);
                if metrics.height > iy {
                    iy = metrics.height;
                }
            }
            img.text(rig + 2, top, f, &FONT_DATA, &txt, tc);
        }

        img.text(lef, top, f, &FONT_DATA, txt, tc);
    }

    ImageBuffer::from_raw(ma, ma, img.take_buffer()).unwrap()
}
