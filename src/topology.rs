use std::ops::Mul;

use crate::element::Element;

mod utils;

#[macro_export]
macro_rules! vec3 {
    ($x:expr, $y:expr, $z:expr) => {
        $crate::topology::Vec3::new($x, $y, $z)
    };
}

#[derive(Clone, Debug, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;

    fn mul(self, rhs: f64) -> Self::Output {
        Vec3 {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

/// A chain within a Topology
#[derive(Clone, Debug, PartialEq)]
#[allow(unused)]
pub(crate) struct Chain {
    /// the index of the chain within its topology
    pub(crate) index: usize,

    /// a user-defined identifier for this Chain
    pub(crate) id: String,

    /// the residues composing the chain
    pub(crate) residues: Vec<Residue>,
}

impl Chain {
    fn new(index: usize, id: String) -> Self {
        Self {
            index,
            id,
            residues: Vec::new(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
#[allow(unused)]
pub(crate) struct Residue {
    /// the name of the residue
    name: String,

    /// the index of the residue within its topology
    index: usize,

    /// the chain this residue belongs to. TODO should surely be a reference, or
    /// at least an index referring to the chain in its owning Topology
    chain: Chain,

    /// a user-defined id for this Residue
    id: String,

    /// a user-defined insertion code for this Residue
    insertion_code: String,

    atoms: Vec<Atom>,
}

impl Residue {
    pub(crate) fn new(
        name: String,
        index: usize,
        chain: Chain,
        id: String,
        insertion_code: String,
    ) -> Self {
        Self {
            name,
            index,
            chain,
            id,
            insertion_code,
            atoms: Vec::new(),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Hash, Eq)]
pub struct Atom {
    pub(crate) name: String,
    pub(crate) element: Element,
    /// the index of the atom within its topology
    pub(crate) index: usize,
    pub(crate) id: usize,
}

impl Atom {
    pub(crate) fn new(
        name: String,
        element: Element,
        index: usize,
        id: usize,
    ) -> Self {
        Self {
            name,
            element,
            index,
            id,
        }
    }
}

/// Topology stores the topological information about a system.
///
/// The structure of a Topology object is similar to that of a PDB file. It
/// consists of a set of Chains (often but not always corresponding to polymer
/// chains). Each Chain contains a set of Residues, and each Residue contains a
/// set of Atoms. In addition, the Topology stores a list of which atom pairs
/// are bonded to each other, and the dimensions of the crystallographic unit
/// cell.
///
/// Atom and residue names should follow the PDB 3.0 nomenclature for all
/// molecules for which one exists.
#[derive(Clone, Debug, PartialEq)]
pub struct Topology {
    pub(crate) chains: Vec<Chain>,
    pub(crate) num_residues: usize,
    pub(crate) num_atoms: usize,
    pub(crate) bonds: Vec<(Atom, Atom)>,
    pub(crate) periodic_box_vectors: Option<()>,
}

impl Topology {
    /// Create a new Topology object
    pub fn new() -> Self {
        Self {
            chains: Vec::new(),
            num_residues: 0,
            num_atoms: 0,
            bonds: Vec::new(),
            periodic_box_vectors: None,
        }
    }

    /// Create a new [Chain] and add it to the [Topology]. Also returns the new
    /// [Chain].
    pub(crate) fn add_chain(&mut self, id: String) -> Chain {
        // TODO id can actually be an option, generate it as
        // str(self.chains.len() + 1)
        let chain = Chain::new(self.chains.len(), id);
        self.chains.push(chain.clone());
        chain
    }

    /// Create a new [Residue], add it to `self`, and return it
    pub(crate) fn add_residue(
        &mut self,
        name: &str,
        mut chain: Chain,
        _to_string: String,
        insertion_code: String,
    ) -> Residue {
        if !chain.residues.is_empty()
            && self.num_residues != chain.residues.last().unwrap().index + 1
        {
            panic!("All residues within a chain must be contiguous");
        }
        let id = (self.num_residues + 1).to_string();
        let residue = Residue::new(
            name.to_owned(),
            self.num_residues,
            chain.clone(),
            id,
            insertion_code,
        );
        self.num_residues += 1;
        chain.residues.push(residue.clone());
        residue
    }

    pub fn atoms(&self) -> impl Iterator<Item = &Atom> {
        self.chains
            .iter()
            .flat_map(|chain| chain.residues.iter())
            .flat_map(|residue| residue.atoms.iter())
    }

    pub fn bonds(&self) -> impl Iterator<Item = &(Atom, Atom)> {
        self.bonds.iter()
    }

    pub(crate) fn add_atom(
        &mut self,
        name: String,
        element: Element,
        r: &mut Residue,
        serial_number: usize,
    ) -> Atom {
        if r.atoms.len() > 0
            && self.num_atoms != r.atoms.last().unwrap().index + 1
        {
            panic!("all atoms within a residue must be contiguous");
        }
        let atom = Atom::new(name, element, self.num_atoms, serial_number);
        self.num_atoms += 1;
        r.atoms.push(atom.clone());
        atom
    }

    pub(crate) fn add_bond(&mut self, bond_1: Atom, bond_2: Atom) {
        self.bonds.push((bond_1, bond_2));
    }
}

impl Default for Topology {
    fn default() -> Self {
        Self::new()
    }
}
