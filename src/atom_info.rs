//! This file is pretty much a periodic table

use ElemGroup::*;

/// Element group on the periodic table
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElemGroup {
    /// Alkali metal
    Alkali,
    /// Alkaline earth metal
    AlkEar,
    /// Transition metal
    Trans,
    /// Post-transition metal
    Poor,
    /// Metalloid
    Metoid,
    /// Nonmetal
    Nonmet,
    /// Halogen
    Halogn,
    /// Noble gas
    Noble,
    /// Rare earth
    RarEar,
}

#[derive(Debug, Clone, Copy)]
pub struct AtomData {
    pub name: &'static str,
    pub sym: &'static str,
    pub mass: f32,
    pub group: ElemGroup
}
impl AtomData {
    pub const fn new(name: &'static str, sym: &'static str, mass: f32, group: ElemGroup) -> Self {
        Self {
            name,
            sym,
            mass,
            group
        }
    }
}

/// All of the data, just index the array by the atomic number
pub static ATOM_DATA: &[AtomData] = &[
    AtomData::new("unknown",       "R",  0.0000, Nonmet),
    AtomData::new("hydrogen",      "H",  1.0078, Nonmet),
    AtomData::new("helium",        "He", 4.0026, Noble ),
    AtomData::new("lithium",       "Li", 6.9410, Alkali),
    AtomData::new("beryllium",     "Be", 9.0122, AlkEar),
    AtomData::new("boron",         "B",  10.811, Metoid),
    AtomData::new("carbon",        "C",  12.011, Nonmet),
    AtomData::new("nitrogen",      "N",  14.007, Nonmet),
    AtomData::new("oxygen",        "O",  15.999, Nonmet),
    AtomData::new("fluorine",      "F",  18.998, Halogn),
    AtomData::new("neon",          "Ne", 20.180, Noble ),
    AtomData::new("sodium",        "Na", 22.990, Alkali),
    AtomData::new("magnesium",     "Mg", 24.305, AlkEar),
    AtomData::new("aluminum",      "Al", 26.982, Poor  ), // yes I'm american, sue me
    AtomData::new("silicon",       "Si", 28.086, Metoid),
    AtomData::new("phosphorus",    "P",  30.974, Nonmet),
    AtomData::new("sulfur",        "S",  32.065, Nonmet),
    AtomData::new("chlorine",      "Cl", 35.453, Halogn),
    AtomData::new("argon",         "Ar", 39.948, Noble ),
    AtomData::new("potassium",     "K",  39.098, Alkali),
    AtomData::new("calcium",       "Ca", 40.078, AlkEar),
    AtomData::new("scandium",      "Sc", 44.956, Trans ),
    AtomData::new("titanium",      "Ti", 47.867, Trans ),
    AtomData::new("vanadium",      "V",  50.942, Trans ),
    AtomData::new("chromium",      "Cr", 51.996, Trans ),
    AtomData::new("manganese",     "Mn", 54.938, Trans ),
    AtomData::new("iron",          "Fe", 55.845, Trans ),
    AtomData::new("cobalt",        "Co", 58.933, Trans ),
    AtomData::new("nickel",        "Ni", 58.693, Trans ),
    AtomData::new("copper",        "Cu", 63.546, Trans ),
    AtomData::new("zinc",          "Zn", 65.380, Trans ),
    AtomData::new("gallium",       "Ga", 69.723, Poor  ),
    AtomData::new("germanium",     "Ge", 72.640, Metoid),
    AtomData::new("arsenic",       "As", 74.992, Metoid),
    AtomData::new("selenium",      "Se", 78.960, Nonmet),
    AtomData::new("bromine",       "Br", 79.904, Halogn),
    AtomData::new("krypton",       "Kr", 83.798, Noble ),
    AtomData::new("rubidium",      "Rb", 85.468, Alkali),
    AtomData::new("strontium",     "Sr", 87.620, AlkEar),
    AtomData::new("yttrium",       "Y",  88.906, Trans ),
    AtomData::new("zirconium",     "Zr", 91.224, Trans ),
    AtomData::new("niobium",       "Nb", 92.906, Trans ),
    AtomData::new("molybdenum",    "Mo", 95.950, Trans ),
    AtomData::new("technetium",    "Tc", 98.000, Trans ),
    AtomData::new("ruthenium",     "Ru", 101.07, Trans ),
    AtomData::new("rhodium",       "Rh", 102.91, Trans ),
    AtomData::new("palladium",     "Pd", 106.42, Trans ),
    AtomData::new("silver",        "Ag", 107.87, Trans ),
    AtomData::new("cadmium",       "Cd", 112.41, Trans ),
    AtomData::new("indium",        "In", 114.82, Poor  ),
    AtomData::new("tin",           "Sn", 118.71, Poor  ),
    AtomData::new("antimony",      "Sb", 121.76, Metoid),
    AtomData::new("tellurium",     "Te", 127.60, Metoid),
    AtomData::new("iodine",        "I",  126.90, Halogn),
    AtomData::new("xenon",         "Xe", 131.29, Noble ),
    AtomData::new("cesium",        "Cs", 132.91, Alkali),
    AtomData::new("barium",        "Ba", 137.33, AlkEar),
    AtomData::new("lanthanum",     "La", 138.91, RarEar),
    AtomData::new("cerium",        "Ce", 140.12, RarEar),
    AtomData::new("praseodymium",  "Pr", 140.91, RarEar),
    AtomData::new("neodymium",     "Nd", 144.24, RarEar),
    AtomData::new("prometheum",    "Pm", 145.00, RarEar),
    AtomData::new("samarium",      "Sm", 150.36, RarEar),
    AtomData::new("europium",      "Eu", 151.96, RarEar),
    AtomData::new("gadolinium",    "Gd", 157.25, RarEar),
    AtomData::new("terbium",       "Tb", 158.93, RarEar),
    AtomData::new("dysprosium",    "Dy", 162.50, RarEar),
    AtomData::new("holmium",       "Ho", 164.93, RarEar),
    AtomData::new("erbium",        "Er", 167.26, RarEar),
    AtomData::new("thulium",       "Tm", 168.93, RarEar),
    AtomData::new("ytterbium",     "Yb", 173.04, RarEar),
    AtomData::new("lutetium",      "Lu", 174.97, RarEar),
    AtomData::new("hafnium",       "Hf", 178.49, Trans ),
    AtomData::new("tantalum",      "Ta", 180.95, Trans ),
    AtomData::new("tungsten",      "W",  183.84, Trans ),
    AtomData::new("rhenium",       "Re", 186.21, Trans ),
    AtomData::new("osmium",        "Os", 190.23, Trans ),
    AtomData::new("iridium",       "Ir", 192.22, Trans ),
    AtomData::new("platinum",      "Pt", 195.08, Trans ),
    AtomData::new("gold",          "Au", 196.97, Trans ),
    AtomData::new("mercury",       "Hg", 200.59, Trans ),
    AtomData::new("thallium",      "Tl", 204.38, Poor  ),
    AtomData::new("lead",          "Pb", 207.20, Poor  ),
    AtomData::new("bismuth",       "Bi", 208.98, Poor  ),
    AtomData::new("polonium",      "Po", 209.00, Poor  ),
    AtomData::new("astatine",      "At", 210.00, Halogn),
    AtomData::new("radon",         "Rn", 222.00, Noble ),
    AtomData::new("francium",      "Fr", 223.00, Alkali),
    AtomData::new("radium",        "Ra", 226.00, AlkEar),
    AtomData::new("actinium",      "Ac", 227.00, RarEar),
    AtomData::new("thorium",       "Th", 232.04, RarEar),
    AtomData::new("protactinium",  "Pa", 231.04, RarEar),
    AtomData::new("uranium",       "U",  238.03, RarEar),
    AtomData::new("neptunium",     "Np", 237.05, RarEar),
    AtomData::new("plutonium",     "Pu", 244.00, RarEar),
    AtomData::new("americium",     "Am", 243.00, RarEar),
    AtomData::new("curium",        "Cm", 247.00, RarEar),
    AtomData::new("berkelium",     "Bk", 247.00, RarEar),
    AtomData::new("californium",   "Cf", 251.00, RarEar),
    AtomData::new("einsteinium",   "Es", 252.00, RarEar),
    AtomData::new("fermium",       "Fm", 257.00, RarEar),
    AtomData::new("mendelvium",    "Md", 258.00, RarEar),
    AtomData::new("nobelium",      "No", 259.00, RarEar),
    AtomData::new("lawrencium",    "Lr", 262.00, RarEar),
    AtomData::new("rutherfordium", "Rf", 267.00, Trans ),
    AtomData::new("dubnium",       "Db", 262.00, Trans ),
    AtomData::new("seaborgium",    "Sg", 269.00, Trans ),
    AtomData::new("bohrium",       "Bh", 264.00, Trans ),
    AtomData::new("hassium",       "Hs", 269.00, Trans ),
    AtomData::new("meitnerium",    "Mt", 278.00, Trans ),
    AtomData::new("darmstadtium",  "Ds", 281.00, Trans ),
    AtomData::new("roentgenium",   "Rg", 282.00, Trans ),
    AtomData::new("copernicium",   "Cn", 285.00, Trans ),
    AtomData::new("nihonium",      "Nh", 286.00, Poor  ),
    AtomData::new("flerovium",     "Fl", 289.00, Poor  ),
    AtomData::new("moscovium",     "Mc", 289.00, Poor  ),
    AtomData::new("livermorium",   "Lv", 293.00, Poor  ),
    AtomData::new("tenessine",     "Ts", 294.00, Halogn),
    AtomData::new("oganesson",     "Og", 294.00, Noble ),
];
