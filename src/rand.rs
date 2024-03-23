use crate::prelude::*;
use gcd::Gcd;
use petgraph::data::{Build, Create};
use rand::prelude::*;
use smallvec::SmallVec;

/// Distribution to create a hydrocarbon's empirical formula
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hydrocarbon(pub usize);
impl Default for Hydrocarbon {
    fn default() -> Self {
        Self(20)
    }
}
impl Distribution<EmpiricalFormula> for Hydrocarbon {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> EmpiricalFormula {
        let c = rng.gen_range(1..self.0);
        let h = rng.gen_range((c / 2)..(c * 2 + 2));
        let mut out = EmpiricalFormula::new();
        out.set_atom(6, c);
        out.set_atom(1, h);
        out
    }
}

#[derive(Clone, Copy)]
enum CationOpt {
    Li,
    Na,
    K,
    Rb,
    Cs,
    Be,
    Mg,
    Ca,
    Sr,
    Ba,
    Al,
    Ti,
    Fe,
    Co,
    Ni,
    Cu,
    Zn,
    Ag,
    Au,
    Nh4,
}
impl CationOpt {
    pub fn charge(self, rng: u8) -> u8 {
        use CationOpt::*;
        match self {
            Li | Na | K | Rb | Cs | Ag | Nh4 => 1,
            Be | Mg | Ca | Sr | Ba | Ni | Zn => 2,
            Al | Au => 3,
            Ti => rng % 3 + 2,
            Fe | Co => rng % 2 + 2,
            Cu => rng % 2 + 1,
        }
    }
    pub fn update_emp(self, emp: &mut EmpiricalFormula) {
        use CationOpt::*;
        let atom = match self {
            Li => 3,
            Na => 11,
            K => 19,
            Rb => 37,
            Cs => 54,
            Be => 4,
            Mg => 12,
            Ca => 20,
            Sr => 38,
            Ba => 56,
            Al => 13,
            Ti => 22,
            Fe => 26,
            Co => 27,
            Ni => 28,
            Cu => 29,
            Zn => 30,
            Ag => 47,
            Au => 79,
            Nh4 => {
                emp.add_atom(1, 4);
                7
            }
        };
        emp.add_atom(atom, 1);
    }
    pub fn nodes_edges(self) -> (usize, usize) {
        (1, 0)
    }
    pub fn update_graph<G: Build<NodeWeight = Atom>>(self, graph: &mut G, rng: u8) {
        use CationOpt::*;
        let mut atom = match self {
            Li => Atom::new(3),
            Na => Atom::new(11),
            K => Atom::new(19),
            Rb => Atom::new(37),
            Cs => Atom::new(54),
            Be => Atom::new(4),
            Mg => Atom::new(12),
            Ca => Atom::new(20),
            Sr => Atom::new(38),
            Ba => Atom::new(56),
            Al => Atom::new(13),
            Ti => Atom::new(22),
            Fe => Atom::new(26),
            Co => Atom::new(27),
            Ni => Atom::new(28),
            Cu => Atom::new(29),
            Zn => Atom::new(30),
            Ag => Atom::new(47),
            Au => Atom::new(79),
            Nh4 => {
                let mut atom = Atom::new(7);
                atom.add_hydrogens(4).unwrap();
                atom
            }
        };
        atom.charge = self.charge(rng) as i8;
        graph.add_node(atom);
    }
}

#[derive(Clone, Copy)]
enum AnionOpt {
    F,
    Cl,
    Br,
    I,
    O,
    Oh,
    No3,
    Co3,
    So3,
    So4,
    HCo3,
    HSo3,
    HSo4,
}
impl AnionOpt {
    pub fn charge(self) -> u8 {
        use AnionOpt::*;
        match self {
            F | Cl | Br | I | No3 | HCo3 | HSo3 | HSo4 | Oh => 1,
            O | Co3 | So3 | So4 => 2,
        }
    }
    pub fn update_emp(self, emp: &mut EmpiricalFormula) {
        use AnionOpt::*;
        let (single, o, h) = match self {
            F => (9, 0, 0),
            Cl => (17, 0, 0),
            Br => (35, 0, 0),
            I => (53, 0, 0),
            O => (8, 0, 0),
            Oh => (8, 0, 1),
            No3 => (7, 3, 0),
            Co3 => (6, 3, 0),
            So3 => (16, 3, 0),
            So4 => (16, 4, 0),
            HCo3 => (6, 3, 1),
            HSo3 => (16, 3, 1),
            HSo4 => (16, 4, 1),
        };
        emp.add_atom(single, 1);
        emp.add_atom(8, o);
        emp.add_atom(1, h);
    }
    pub fn nodes_edges(self) -> (usize, usize) {
        use AnionOpt::*;
        match self {
            F | Cl | Br | I | O | Oh => (1, 0),
            No3 | Co3 | So3 | HCo3 | HSo3 => (4, 3),
            So4 | HSo4 => (5, 4),
        }
    }
    pub fn update_graph<G: Build<NodeWeight = Atom, EdgeWeight = Bond>>(self, graph: &mut G) {
        use AnionOpt::*;
        let om = Atom {
            charge: -1,
            ..Atom::new(8)
        };
        let mut oh = Atom::new(8);
        oh.add_hydrogens(1).unwrap();
        match self {
            F | Cl | Br | I => {
                graph.add_node(Atom {
                    charge: -1,
                    ..Atom::new(match self {
                        F => 9,
                        Cl => 17,
                        Br => 35,
                        I => 53,
                        _ => unreachable!(),
                    })
                });
            }
            O => {
                graph.add_node(om);
            }
            Oh => {
                graph.add_node(oh);
            }
            No3 => {
                let c = graph.add_node(Atom {
                    charge: 1,
                    ..Atom::new(7)
                });
                let mut n = graph.add_node(om);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(om);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(Atom::new(8));
                graph.add_edge(c, n, Bond::Double).unwrap();
            }
            Co3 => {
                let c = graph.add_node(Atom::new(6));
                let mut n = graph.add_node(om);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(om);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(Atom::new(8));
                graph.add_edge(c, n, Bond::Double).unwrap();
            }
            So3 => {
                let c = graph.add_node(Atom::new(16));
                let mut n = graph.add_node(om);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(om);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(Atom::new(8));
                graph.add_edge(c, n, Bond::Double).unwrap();
            }
            So4 => {
                let c = graph.add_node(Atom::new(16));
                let mut n = graph.add_node(om);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(om);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(Atom::new(8));
                graph.add_edge(c, n, Bond::Double).unwrap();
                n = graph.add_node(Atom::new(8));
                graph.add_edge(c, n, Bond::Double).unwrap();
            }
            HCo3 => {
                let c = graph.add_node(Atom::new(6));
                let mut n = graph.add_node(oh);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(om);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(Atom::new(8));
                graph.add_edge(c, n, Bond::Double).unwrap();
            }
            HSo3 => {
                let c = graph.add_node(Atom::new(16));
                let mut n = graph.add_node(oh);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(om);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(Atom::new(8));
                graph.add_edge(c, n, Bond::Double).unwrap();
            }
            HSo4 => {
                let c = graph.add_node(Atom::new(16));
                let mut n = graph.add_node(oh);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(om);
                graph.add_edge(c, n, Bond::Single).unwrap();
                n = graph.add_node(Atom::new(8));
                graph.add_edge(c, n, Bond::Double).unwrap();
                n = graph.add_node(Atom::new(8));
                graph.add_edge(c, n, Bond::Double).unwrap();
            }
        };
    }
}

/// Distribution to create a metal salt. Can create any molecule-like graph or an empirical formula.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Salt {
    /// Whether alkali metal salts should be generated
    pub alkali: bool,
    /// Wjhether alkaline earth metal salts should be generated
    pub alkear: bool,
    /// Whether transition metal salts should be generated
    pub trans: bool,
    /// Whether polyatomic cations should be generated
    pub poly_cation: bool,
    /// Whether halids should be generated
    pub halide: bool,
    /// Whether oxides should be generated
    pub oxide: bool,
    /// Whether polyatomic anions should be generated
    pub poly_anion: bool,
}
impl Salt {
    fn select_cation<R: Rng + ?Sized>(self, rng: &mut R) -> CationOpt {
        use CationOpt::*;
        let mut opts = SmallVec::<CationOpt, 18>::new();
        if self.alkali {
            opts.extend_from_slice(&[Li, Na, K, Rb, Cs]);
        }
        if self.alkear {
            opts.extend_from_slice(&[Be, Mg, Ca, Sr, Ba]);
        }
        if self.trans {
            opts.extend_from_slice(&[Al, Ti, Fe, Co, Ni, Cu, Zn, Ag, Au]);
        }
        if self.poly_cation {
            opts.extend_from_slice(&[Nh4]);
        }
        if opts.is_empty() {
            panic!("no cations are available");
        }
        *opts.choose(rng).unwrap()
    }
    fn select_anion<R: Rng + ?Sized>(self, rng: &mut R) -> AnionOpt {
        use AnionOpt::*;
        let mut opts = SmallVec::<AnionOpt, 12>::new();
        if self.halide {
            opts.extend_from_slice(&[F, Cl, Br, I]);
        }
        if self.oxide {
            opts.push(O);
        }
        if self.poly_anion {
            opts.extend_from_slice(&[Oh, No3, Co3, So3, So4, HCo3, HSo3, HSo4]);
        }
        if opts.is_empty() {
            panic!("no anions are available");
        }
        *opts.choose(rng).unwrap()
    }
}
impl Default for Salt {
    fn default() -> Self {
        Self {
            alkali: true,
            alkear: true,
            trans: true,
            oxide: true,
            halide: true,
            poly_cation: true,
            poly_anion: true,
        }
    }
}
impl Distribution<EmpiricalFormula> for Salt {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> EmpiricalFormula {
        let cation = self.select_cation(rng);
        let anion = self.select_anion(rng);
        let r = rng.gen();
        let cc = cation.charge(r);
        let ac = anion.charge();
        let gcd = cc.gcd(ac);
        let mut out = EmpiricalFormula::new();
        for _ in 0..(ac / gcd) {
            cation.update_emp(&mut out);
        }
        for _ in 0..(cc / gcd) {
            anion.update_emp(&mut out);
        }
        out
    }
}
impl<G: Create<NodeWeight = Atom, EdgeWeight = Bond>> Distribution<G> for Salt {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> G {
        let cation = self.select_cation(rng);
        let anion = self.select_anion(rng);
        let r = rng.gen();
        let cc = cation.charge(r);
        let ac = anion.charge();
        let gcd = cc.gcd(ac);
        let nc = ac / gcd;
        let na = cc / gcd;
        let (cn, ce) = cation.nodes_edges();
        let (an, ae) = anion.nodes_edges();
        let mut out = G::with_capacity(
            cn * nc as usize + an * na as usize,
            ce * nc as usize + ae * na as usize,
        );
        for _ in 0..(ac / gcd) {
            cation.update_graph(&mut out, r);
        }
        for _ in 0..(cc / gcd) {
            anion.update_graph(&mut out);
        }
        out
    }
}
