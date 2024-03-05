use super::arena::*;
use super::*;
pub mod graph_traits;
mod node_impls;
pub use node_impls::*;

/// A `Molecule` acts like a graph, and can have graph algorithms used on it. It's immutable, with all
/// mutations making (efficient) copies.
#[derive(Debug, Clone, Copy)]
pub struct Molecule<Ix, R> {
    arena: R,
    index: Ix,
}
impl<Ix, R> Molecule<Ix, R> {
    pub fn from_arena<'a, 'b: 'a, A: ArenaAccessible<Ix = Ix, Access<'a> = R> + 'a>(
        arena: &'b A,
        index: Ix,
    ) -> Self
    where
        Self: 'a,
    {
        Self {
            arena: arena.get_accessor(),
            index,
        }
    }

    pub fn from_mut_arena<'a, 'b: 'a, A: ArenaAccessibleMut<Ix = Ix, AccessMut<'a> = R> + 'a>(
        arena: &'b A,
        index: Ix,
    ) -> Self
    where
        Self: 'a,
    {
        Self {
            arena: arena.get_accessor_mut(),
            index,
        }
    }

    pub fn arena(&self) -> R::Ref<'_>
    where
        R: ArenaAccessor<Ix = Ix>,
    {
        self.arena.get_arena()
    }

    pub fn arena_mut(&self) -> R::RefMut<'_>
    where
        R: ArenaAccessorMut<Ix = Ix>,
    {
        self.arena.get_arena_mut()
    }
}

impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> Molecule<Ix, R> {
    /// Check if this molecule contains the underlying group.
    pub fn contains(&self, group: Ix) -> bool {
        self.arena.get_arena().contains_group(self.index, group)
    }

    /// Get an atom in this molecule. Returns a value because of possible `MolRepr::Modify`s.
    pub fn get_atom(&self, mut idx: NodeIndex<Ix>) -> R::MappedRef<'_, Atom> {
        R::map_ref(self.arena.get_arena(), |arena| {
            let mut ix = self.index;
            loop {
                match &arena.parts[ix.index()] {
                    (MolRepr::Redirect(i), _) => ix = *i,
                    (MolRepr::Modify(m), _) => {
                        if let Some(a) = m.patch.get(&idx.0) {
                            return a;
                        } else {
                            ix = m.base
                        }
                    }
                    _ => todo!(),
                }
            }
        })
    }
}

impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>, S: ArenaAccessor<Ix = Ix>> PartialEq<Molecule<Ix, S>>
    for Molecule<Ix, R>
{
    fn eq(&self, other: &Molecule<Ix, S>) -> bool {
        self.index == other.index && std::ptr::eq(&*self.arena(), &*other.arena())
    }
}
