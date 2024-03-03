use crate::atom_info::*;
use crate::core::*;
use petgraph::dot::*;
use petgraph::visit::*;
use petgraph::Undirected;
use std::fmt::{Debug, Display};

pub fn dot_atom_color(num: u8) -> &'static str {
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

pub fn bond_style<
    G: GraphRef + Data<NodeWeight = Atom, EdgeWeight = Bond>,
    E: EdgeRef<Weight = Bond>,
>(
    _: G,
    bond: E,
) -> String {
    match *bond.weight() {
        Bond::Single => r#"color="black""#,
        Bond::Double => r#"color="black:black""#,
        Bond::Triple => r#"color="black:black:black""#,
        Bond::Quad => r#"color="black:black:black:black""#,
        Bond::Aromatic => r#"color="black:black", style="dashed""#,
        Bond::Left => r#"color="black", label="left""#,
        Bond::Right => r#"color="black", label="right""#,
        Bond::Non => r#"color="black", style="invis""#,
        _ => panic!("invalid bond!"),
    }
    .to_string()
}
pub fn atom_style<G: GraphRef + Data<NodeWeight = Atom, EdgeWeight = Bond> + IntoNodeReferences>(
    _: G,
    atom: G::NodeRef,
) -> String {
    let num = atom.weight().protons;
    format!("color={:?}", dot_atom_color(num))
}

pub fn fmt_as_dot<'a, G>(graph: G) -> impl Debug + Display + 'a
where
    G: GraphRef
        + GraphProp<EdgeType = Undirected>
        + Data<NodeWeight = Atom, EdgeWeight = Bond>
        + IntoNodeReferences
        + IntoEdgeReferences
        + NodeIndexable
        + 'a,
{
    Dot::with_attr_getters(graph, &[Config::EdgeNoLabel], &bond_style, &atom_style)
}
