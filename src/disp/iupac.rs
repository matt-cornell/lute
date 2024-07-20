bitflags::bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub struct IupacConfig: u8 {
        /// No substitutions, just the PIN.
        const PREFERRED = 0;
        /// Use iso-, neo-, and tert- prefixes
        const ISO_NEO_TERT = 1;
        /// Don't generate names like "dimethyl ether", prefer "methoxymethane"
        const NEVER_ETHER = 2;
        /// Give chemicals more common names if possible
        const COMMON_SUBS = 4;
        /// Give names like "formic acid", forces -aldehyde suffix
        const FORM_ACET = 8;
        /// Combination of options to give the most common names
        const COMMON_NAME = 13;
    }
}

pub fn iupac_name<G>(graph: G, cfg: IupacConfig) -> String {
    "TODO".to_string()
}
