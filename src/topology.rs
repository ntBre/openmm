#[macro_export]
macro_rules! vec3 {
    ($x:expr, $y:expr, $z:expr) => {
        $crate::Vec3::new($x, $y, $z)
    };
}

#[derive(Clone)]
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

/// A chain within a Topology
#[derive(Clone)]
#[allow(unused)]
pub(crate) struct Chain {
    /// the index of the chain within its topology
    index: usize,

    /// the topology this chain belongs to. TODO surely this needs to be a
    /// reference, but that's very scary
    topology: Topology,

    /// a user-defined identifier for this Chain
    id: String,

    /// the residues composing the chain
    residues: Vec<Residue>,
}

impl Chain {
    fn new(index: usize, topology: Topology, id: String) -> Self {
        Self {
            index,
            topology,
            id,
            residues: Vec::new(),
        }
    }
}

#[derive(Clone)]
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

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Atom {
    pub(crate) index: usize,
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
#[derive(Clone)]
#[allow(unused)]
pub struct Topology {
    chains: Vec<Chain>,
    num_residues: usize,
    num_atoms: usize,
    bonds: Vec<(Atom, Atom)>,
    periodic_box_vectors: Option<()>,
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
        let chain = Chain::new(self.chains.len(), self.clone(), id);
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

    pub(crate) fn bonds(&self) -> impl Iterator<Item = &(Atom, Atom)> {
        self.bonds.iter()
    }
}

impl Default for Topology {
    fn default() -> Self {
        Self::new()
    }
}
