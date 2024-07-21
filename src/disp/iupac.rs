use crate::core::*;
use crate::graph::misc::DataValueMap;
use crate::graph::{BitFiltered, ConnectedGraphIter};
use petgraph::visit::*;
use smallvec::SmallVec;

bitflags::bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub struct IupacConfig: u8 {
        /// No substitutions, just the PIN.
        const PREFERRED = 0;
        /// Use iso-, neo-, and tert- prefixes
        const ISO_NEO_TERT = 1;
        /// Don't generate names like "dimethyl ether", prefer "methoxymethane".
        const NEVER_ETHER = 2;
        /// Give chemicals more common names if possible.
        const COMMON_SUBS = 4;
        /// Give names like "formic acid", forces -aldehyde suffix for aldehydes.
        const FORM_ACET = 8;
        /// Don't have distinct R/S and E/Z stereoisomers.
        const NO_STEREO = 16;
        /// No numbers, e.g. 1-propanol and 2-propanol both become "propanol".
        /// "Isopropanol" could still be generated with this set.
        const NO_NUMBERS = 32;
        /// Generate "hyrochloric acid" for HCl instead of "hydrogen chloride"
        const HYDRO_ACIDS = 64;
        /// Combination of options to give the most common names.
        const COMMON_NAME = 29;
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
enum Halogen {
    Hal,
    At,
    Br,
    Cl,
    F,
    I,
}
impl Halogen {
    pub fn from_u8(p: u8) -> Self {
        match p {
            0 => Self::Hal,
            9 => Self::F,
            17 => Self::Cl,
            35 => Self::Br,
            53 => Self::I,
            85 => Self::At,
            _ => unreachable!(),
        }
    }
    pub fn fragment(self) -> &'static str {
        match self {
            Self::Hal => "hal",
            Self::At => "astat",
            Self::Br => "brom",
            Self::Cl => "chlor",
            Self::F => "fluor",
            Self::I => "iod",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Feature<N> {
    Alcohol {
        o: N,
        base: bool,
    },
    Ether(N),
    Ketone {
        c: N,
        o: N,
    },
    Thiol {
        s: N,
        base: bool,
    },
    Thioether(N),
    Thioketone {
        c: N,
        s: N,
    },
    CarbAcid {
        c: N,
        os: [N; 2],
        base: bool,
    },
    Ester {
        c: N,
        eth_o: N,
        ket_o: N,
    },
    Thioester {
        c: N,
        s: N,
        o: N,
    },
    Amine(N, i8),
    Amide {
        n: N,
        c: N,
        o: N,
    },
    Imine {
        c: N,
        n: N,
    },
    Imide {
        n: N,
        c1: N,
        o1: N,
        c2: N,
        o2: N,
    },
    Halide {
        hal: Halogen,
        c: N,
        x: N,
    },
    HypoHalide {
        hal: Halogen,
        c: N,
        o: N,
        x: N,
    },
    DioxyHalide {
        hal: Halogen,
        c: N,
        o: N,
        x: N,
        os: N,
    },
    TrioxyHalide {
        hal: Halogen,
        c: N,
        o: N,
        x: N,
        os: [N; 2],
    },
    PeroxyHalide {
        hal: Halogen,
        c: N,
        o: N,
        x: N,
        os: [N; 3],
    },
    AcylHalide {
        hal: Halogen,
        c: N,
        x: N,
        o: N,
    },
    Sulfate {
        s: N,
        o: N,
        os: [N; 2],
    },
    Phenyl([N; 6]),
    Benzyl(N, [N; 6]),
    Benzoyl {
        c: N,
        o: N,
        ring: [N; 6],
    },
}

#[derive(Debug, Clone, Copy)]
enum SubKind {
    Mol,
    Base,
    R,
    H,
}

#[derive(Debug, Clone, Copy)]
pub enum InvalidMolecule<N> {
    InvalidValence {
        atom: N,
        num: isize,
        expected: &'static [u8],
    },
    InvalidCharge {
        atom: N,
        charge: i8,
        expected: &'static [i8],
    },
}

pub fn iupac_name<
    G: DataValueMap<NodeWeight = Atom, EdgeWeight = Bond>
        + IntoNodeReferences
        + NodeIndexable
        + IntoEdges,
>(
    graph: G,
    cfg: IupacConfig,
) -> Result<String, InvalidMolecule<G::NodeId>> {
    let mut features = Vec::new();
    let mut frags: SmallVec<_, 3> = ConnectedGraphIter::new(graph)
        .iter(graph)
        .map(|filter| {
            features.clear();
            let graph = BitFiltered::<G, usize, 1, false>::new(graph, filter);
            for n in graph.node_references() {
                let atom = *n.weight();
                match atom.protons {
                    7 => {}
                    8 => {
                        if atom.data.other() == 0 {
                            let total =
                                (atom.data.hydrogen() + atom.data.single() + atom.data.unknown())
                                    as i8
                                    - atom.charge;
                            if total != 2 {
                                return Err(InvalidMolecule::InvalidValence {
                                    atom: n.id(),
                                    num: total as isize,
                                    expected: &[2],
                                });
                            }
                            match atom.data.hydrogen() {
                                0 => match atom.charge {
                                    0 => {
                                        if atom.data.unknown() == 2 {
                                            return Ok(("ether".to_string(), 0));
                                        } else {
                                            features.push(Feature::Ether(n.id()));
                                        }
                                    }
                                    -1 => {
                                        if atom.data.unknown() == 1 {
                                            return Ok(("oxide".to_string(), -1));
                                        } else {
                                            features.push(Feature::Alcohol {
                                                o: n.id(),
                                                base: true,
                                            })
                                        }
                                    }
                                    _ => {}
                                },
                                1 => match atom.charge {
                                    0 => {
                                        if atom.data.unknown() == 1 {
                                            return Ok(("hydroxyl".to_string(), -1));
                                        } else {
                                            features.push(Feature::Alcohol {
                                                o: n.id(),
                                                base: false,
                                            })
                                        }
                                    }
                                    _ => {}
                                },
                                2 => {
                                    if atom.charge == 0 {
                                        return Ok(("water".to_string(), 0));
                                    }
                                }
                                3 => {
                                    if atom.charge == 1 {
                                        return Ok(("hydronium".to_string(), 1));
                                    }
                                }
                                _ => {}
                            }
                        } else if atom.data.other() == 1
                            && atom.data.unknown() == 0
                            && atom.data.single() == 0
                            && atom.data.hydrogen() == 0
                        {
                            let Some(ed) = graph.edges(n.id()).next() else {
                                continue;
                            };
                            let neighbor = graph.node_weight(ed.target()).unwrap();
                            if neighbor.protons == 6 {
                                if neighbor.data.other() == 1
                                    && neighbor.data.unknown() == 0
                                    && neighbor.data.single() == 0
                                    && neighbor.data.hydrogen() == 0
                                {
                                    return Ok(("carbon monoxide".to_string(), 0));
                                } else if *ed.weight() == Bond::Double {
                                    features.push(Feature::Ketone {
                                        c: ed.target(),
                                        o: n.id(),
                                    });
                                }
                            }
                        }
                    }
                    16 => {
                        if atom.data.other() == 0 {
                            let total =
                                (atom.data.hydrogen() + atom.data.single() + atom.data.unknown())
                                    as i8
                                    - atom.charge;
                            if total != 2 {
                                return Err(InvalidMolecule::InvalidValence {
                                    atom: n.id(),
                                    num: total as isize,
                                    expected: &[2],
                                });
                            }
                            match atom.data.hydrogen() {
                                0 => match atom.charge {
                                    0 => {
                                        if atom.data.unknown() == 2 {
                                            return Ok(("thioether".to_string(), 0));
                                        } else {
                                            features.push(Feature::Thioether(n.id()));
                                        }
                                    }
                                    -1 => {
                                        if atom.data.unknown() == 1 {
                                            return Ok(("sulfide".to_string(), -1));
                                        } else {
                                            features.push(Feature::Thiol {
                                                s: n.id(),
                                                base: true,
                                            })
                                        }
                                    }
                                    _ => {}
                                },
                                1 => match atom.charge {
                                    0 => {
                                        if atom.data.unknown() == 1 {
                                            return Ok(("thiol".to_string(), -1));
                                        } else {
                                            features.push(Feature::Thiol {
                                                s: n.id(),
                                                base: false,
                                            })
                                        }
                                    }
                                    _ => {}
                                },
                                2 => {
                                    if atom.charge == 0 {
                                        return Ok((
                                            if cfg.contains(IupacConfig::HYDRO_ACIDS) {
                                                "hydrosulfuric acid"
                                            } else {
                                                "hydrogen sulfide"
                                            }
                                            .to_string(),
                                            0,
                                        ));
                                    }
                                }
                                3 => {
                                    if atom.charge == 1 {
                                        return Ok(("sulfonium".to_string(), 1));
                                    }
                                }
                                _ => {}
                            }
                        } else if atom.data.other() == 1
                            && atom.data.unknown() == 0
                            && atom.data.single() == 0
                            && atom.data.hydrogen() == 0
                        {
                            let Some(ed) = graph.edges(n.id()).next() else {
                                continue;
                            };
                            let neighbor = graph.node_weight(ed.target()).unwrap();
                            if neighbor.protons == 6 {
                                if neighbor.data.other() == 1
                                    && neighbor.data.unknown() == 0
                                    && neighbor.data.single() == 0
                                    && neighbor.data.hydrogen() == 0
                                {
                                    return Ok(("sulfur monoxide".to_string(), 0));
                                } else if *ed.weight() == Bond::Double {
                                    features.push(Feature::Thioketone {
                                        c: ed.target(),
                                        s: n.id(),
                                    });
                                }
                            }
                        }
                    }
                    p @ (0 | 9 | 17 | 35 | 53 | 85) if p != 0 || atom.isotope == 0x0100 => {
                        let hal = Halogen::from_u8(p);
                        let frag = hal.fragment();
                        if atom.data.other() == 0 {
                            match (
                                atom.data.single(),
                                atom.data.hydrogen(),
                                atom.data.unknown(),
                            ) {
                                (1, 0, 0) => {
                                    let n1 = graph.neighbors(n.id()).next().unwrap();
                                    let a = graph.node_weight(n1).unwrap();
                                    if a.protons == 8 {
                                        if a.data.other() != 0 {
                                            let valence = graph
                                                .edges(n1)
                                                .map(|i| i.weight().bond_count())
                                                .sum::<f32>()
                                                as isize
                                                - a.charge as isize;
                                            if valence == 2 {
                                                return Err(InvalidMolecule::InvalidCharge {
                                                    atom: n1,
                                                    charge: a.charge,
                                                    expected: &[0, -1],
                                                });
                                            } else {
                                                return Err(InvalidMolecule::InvalidValence {
                                                    atom: n1,
                                                    num: valence,
                                                    expected: &[2],
                                                });
                                            }
                                        }
                                        match (a.data.single(), a.data.hydrogen(), a.data.unknown())
                                        {
                                            (2, 0, 0) => {
                                                let n2 = graph
                                                    .neighbors(n1)
                                                    .find(|&n2| n2 != n.id())
                                                    .unwrap();
                                                features.push(Feature::HypoHalide {
                                                    hal,
                                                    c: n2,
                                                    o: n1,
                                                    x: n.id(),
                                                });
                                            }
                                            (1, 1, 0) => {
                                                return Ok((format!("hypo{frag}ous acid"), 0));
                                            }
                                            (1, 0, u @ (0 | 1)) => {
                                                let valence = u as i8 - atom.charge;
                                                if valence == 1 {
                                                    return Ok((
                                                        format!("hypo{frag}ite"),
                                                        atom.charge as _,
                                                    ));
                                                } else {
                                                    return Err(InvalidMolecule::InvalidValence {
                                                        atom: n.id(),
                                                        num: valence as _,
                                                        expected: &[1],
                                                    });
                                                }
                                            }
                                            (a, b, c) => {
                                                let valence = (a + b + c) as i8 - atom.charge;
                                                if valence == 2 {
                                                    return Err(InvalidMolecule::InvalidCharge {
                                                        atom: n.id(),
                                                        charge: atom.charge,
                                                        expected: &[0, -1],
                                                    });
                                                } else {
                                                    return Err(InvalidMolecule::InvalidValence {
                                                        atom: n.id(),
                                                        num: valence as _,
                                                        expected: &[1],
                                                    });
                                                }
                                            }
                                        }
                                    } else {
                                        features.push(Feature::Halide {
                                            hal,
                                            c: n1,
                                            x: n.id(),
                                        });
                                    }
                                    continue;
                                }
                                (0, 1, 0) => {
                                    if cfg.contains(IupacConfig::HYDRO_ACIDS) {
                                        return Ok((format!("hyrdro{frag}ic acid"), 0));
                                    } else {
                                        return Ok((format!("hydrogen {frag}ide"), 0));
                                    }
                                }
                                (0, 0, u @ (0 | 1)) => {
                                    let valence = u as i8 - atom.charge;
                                    if valence == 1 {
                                        return Ok((format!("{frag}ide"), atom.charge as _));
                                    } else {
                                        return Err(InvalidMolecule::InvalidValence {
                                            atom: n.id(),
                                            num: valence as _,
                                            expected: &[1],
                                        });
                                    }
                                }
                                (a, b, c) => {
                                    let valence = (a + b + c) as i8 - atom.charge;
                                    if valence == 2 {
                                        return Err(InvalidMolecule::InvalidCharge {
                                            atom: n.id(),
                                            charge: atom.charge,
                                            expected: &[0, -1],
                                        });
                                    } else {
                                        return Err(InvalidMolecule::InvalidValence {
                                            atom: n.id(),
                                            num: valence as _,
                                            expected: &[1],
                                        });
                                    }
                                }
                            }
                        } else if atom.data.single() == 1
                            && atom.data.hydrogen() == 0
                            && atom.data.unknown() == 0
                            && atom.data.other() <= 3
                        {
                            let check = 'check: {
                                let mut kind = None;
                                for edge in graph.edges(n.id()) {
                                    let other = graph.node_weight(edge.target()).unwrap();
                                    if other.protons != 8 {
                                        break 'check None;
                                    }
                                    match *edge.weight() {
                                        Bond::Double => {}
                                        Bond::Single if kind.is_none() => {
                                            let o = edge.target();
                                            match (
                                                other.data.hydrogen(),
                                                other.data.single(),
                                                other.data.unknown(),
                                                other.charge,
                                            ) {
                                                (1, 1, 0, 0) => kind = Some((SubKind::H, o)),
                                                (0, 2, 0, 0) => kind = Some((SubKind::Mol, o)),
                                                (0, 1, 1, 0) => kind = Some((SubKind::R, o)),
                                                (0, 1, 0, -1) => kind = Some((SubKind::Base, o)),
                                                _ => {}
                                            }
                                        }
                                        _ => break 'check None,
                                    }
                                }
                                kind
                            };
                            if let Some((kind, o)) = check {
                                let neighbor =
                                    || graph.neighbors(o).find(|&i| i != n.id()).unwrap();
                                let mut os = graph.edges(n.id()).filter_map(|e| {
                                    (*e.weight() == Bond::Double).then_some(e.target())
                                });
                                match (atom.data.other(), kind) {
                                    (1, SubKind::H) => return Ok((format!("{frag}ous acid"), 0)),
                                    (1, SubKind::R) => return Ok((format!("{frag}ite"), 0)),
                                    (1, SubKind::Base) => return Ok((format!("{frag}ite"), -1)),
                                    (1, SubKind::Mol) => {
                                        features.push(Feature::DioxyHalide {
                                            hal,
                                            c: neighbor(),
                                            o,
                                            x: n.id(),
                                            os: os.next().unwrap(),
                                        });
                                    }
                                    (2, SubKind::H) => return Ok((format!("{frag}ic acid"), 0)),
                                    (2, SubKind::R) => return Ok((format!("{frag}ate"), 0)),
                                    (2, SubKind::Base) => return Ok((format!("{frag}ate"), -1)),
                                    (2, SubKind::Mol) => {
                                        features.push(Feature::TrioxyHalide {
                                            hal,
                                            c: neighbor(),
                                            o,
                                            x: n.id(),
                                            os: os.collect::<Vec<_>>().try_into().ok().unwrap(),
                                        });
                                    }
                                    (3, SubKind::H) => return Ok((format!("per{frag}ic acid"), 0)),
                                    (3, SubKind::R) => return Ok((format!("per{frag}ate"), 0)),
                                    (3, SubKind::Base) => return Ok((format!("per{frag}ate"), -1)),
                                    (3, SubKind::Mol) => {
                                        features.push(Feature::PeroxyHalide {
                                            hal,
                                            c: neighbor(),
                                            o,
                                            x: n.id(),
                                            os: os.collect::<Vec<_>>().try_into().ok().unwrap(),
                                        });
                                    }
                                    _ => unreachable!(),
                                }
                            }
                        }
                    }
                    _ => {}
                }
            }
            todo!()
        })
        .collect::<Result<_, _>>()?;
    frags.sort_by_key(|x| i32::MAX - x.1);
    if let Some(((out, _), rest)) = frags.split_first_mut() {
        for chunk in rest.chunks(8) {
            out.reserve(chunk.iter().map(|i| i.0.len()).sum::<usize>() + chunk.len());
            for (n, _) in chunk {
                out.push(' ');
                out.push_str(n);
            }
        }
        Ok(std::mem::take(out))
    } else {
        Ok(String::new())
    }
}
