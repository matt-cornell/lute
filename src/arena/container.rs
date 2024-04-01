use ahash::HashMap;

/// A `Container` corresponds to a reaction vessel. It's at this layer that actual reactions are
/// handled.
#[derive(Debug, Clone)]
pub struct Container<Ix> {
    pub quantities: HashMap<Ix, f64>,
    pub temperature: f64,
}
