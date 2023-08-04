#![feature(lazy_cell)]
#![allow(unused)]

use std::{
    collections::HashMap,
    error::Error,
    io::{BufRead, BufReader, Read},
    path::Path,
};

use element::{Element, BY_SYMBOL, EP};

macro_rules! q {
    ($v:expr, $u:expr) => {
        Quantity::new($v, $u)
    };
}

mod element;

#[derive(Clone)]
pub struct Quantity<T> {
    pub value: T,
    pub unit: String,
}

impl<T> Quantity<T> {
    pub fn new(value: T, unit: impl Into<String>) -> Self {
        Self {
            value,
            unit: unit.into(),
        }
    }
}

#[allow(unused)]
macro_rules! vec3 {
    ($x:expr, $y:expr, $z:expr) => {
        $crate::Vec3::new($x, $y, $z)
    };
}

pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

struct Location {
    alternate_location_indicator: String,
    xyz: Quantity<Vec3>,
    occupancy: f64,
    temperature_factor: Quantity<f64>,
    residue_name_with_spaces: String,
}

impl Location {
    fn new(
        alternate_location_indicator: String,
        xyz: Quantity<Vec3>,
        occupancy: f64,
        temperature_factor: Quantity<f64>,
        residue_name_with_spaces: String,
    ) -> Self {
        Self {
            alternate_location_indicator,
            xyz,
            occupancy,
            temperature_factor,
            residue_name_with_spaces,
        }
    }
}

/// Atom represents one atom in a PDB structure.
#[derive(Default)]
struct Atom {
    is_first_atom_in_chain: bool,
    is_final_atom_in_chain: bool,
    is_first_residue_in_chain: bool,
    is_final_residue_in_chain: bool,
    record_name: String,
    serial_number: usize,
    name_with_spaces: String,
    residue_name_with_spaces: String,
    residue_name: String,
    chain_id: String,
    residue_number: usize,
    insertion_code: String,
    locations: HashMap<String, Location>,
    default_location_id: String,
    segment_id: String,
    element_symbol: String,
    formal_charge: Option<isize>,
    element: Option<Element>,
    model_number: usize,
}

impl Atom {
    /// Create a new pdb.Atom from an ATOM or HETATM line.
    ///
    /// Example line:
    /// ATOM   2209  CB  TYR A 299       6.167  22.607  20.046  1.00  8.12           C
    /// 00000000011111111112222222222333333333344444444445555555555666666666677777777778
    /// 12345678901234567890123456789012345678901234567890123456789012345678901234567890
    ///
    /// ATOM line format description from
    ///   http://deposit.rcsb.org/adit/docs/pdb_atom_format.html:
    ///
    /// COLUMNS        DATA TYPE       CONTENTS
    /// --------------------------------------------------------------------------------
    ///  1 -  6        Record name     "ATOM  "
    ///  7 - 11        Integer         Atom serial number.
    /// 13 - 16        Atom            Atom name.
    /// 17             Character       Alternate location indicator.
    /// 18 - 20        Residue name    Residue name.
    /// 22             Character       Chain identifier.
    /// 23 - 26        Integer         Residue sequence number.
    /// 27             AChar           Code for insertion of residues.
    /// 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
    /// 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
    /// 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
    /// 55 - 60        Real(6.2)       Occupancy (Default = 1.0).
    /// 61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
    /// 73 - 76        LString(4)      Segment identifier, left-justified.
    /// 77 - 78        LString(2)      Element symbol, right-justified.
    /// 79 - 80        LString(2)      Charge on the atom.
    fn new(pdb_line: &str, pdbstructure: &mut PdbStructure) -> Self {
        let mut ret = Self::default();
        ret.record_name = pdb_line[..6].trim().to_owned();

        // TODO this could also be in hex supposedly
        ret.serial_number = match pdb_line[6..11].parse() {
            Ok(i) => i,
            Err(_) => pdbstructure.next_atom_number,
        };

        ret.name_with_spaces = pdb_line[12..16].to_owned();
        let alternate_location_indicator = &pdb_line[16..17];

        ret.residue_name_with_spaces = pdb_line[17..20].to_owned();

        // in some MD codes, notably ffamber in gromacs, residue name has a
        // fourth character in column 21
        let possible_fourth_character = &pdb_line[20..21];
        if possible_fourth_character != " " {
            if ret.residue_name_with_spaces.trim().len() != 3 {
                panic!("misaligned residue name: {pdb_line}");
            }
            ret.residue_name_with_spaces += possible_fourth_character;
        }
        ret.residue_name = ret.residue_name_with_spaces.trim().to_owned();

        ret.chain_id = pdb_line[21..22].to_owned();

        ret.residue_number = pdb_line[22..26]
            .parse()
            .expect("failed to parse residue number");

        ret.insertion_code = pdb_line[26..27].to_owned();

        let x: f64 = pdb_line[30..38].parse().expect("failed to parse x coord");
        let y: f64 = pdb_line[38..46].parse().expect("failed to parse y coord");
        let z: f64 = pdb_line[46..54].parse().expect("failed to parse z coord");

        let occupancy = pdb_line[54..60].parse().unwrap_or(1.0);

        let temperature_factor = match pdb_line[60..66].parse::<f64>() {
            Ok(t) => Quantity::new(t, "angstroms**2"),
            Err(_) => Quantity::new(0.0, "angstroms**2"),
        };

        let loc = Location::new(
            alternate_location_indicator.to_owned(),
            Quantity::new(Vec3::new(x, y, z), "angstroms"),
            occupancy,
            temperature_factor,
            ret.residue_name_with_spaces.clone(),
        );

        ret.locations
            .insert(alternate_location_indicator.to_owned(), loc);
        ret.default_location_id = alternate_location_indicator.to_owned();

        ret.segment_id = pdb_line[72..76].trim().to_owned();
        ret.element_symbol = pdb_line[76..78].trim().to_owned();

        ret.formal_charge = match pdb_line[78..80].parse() {
            Ok(c) => Some(c),
            Err(_) => None,
        };

        if ret.element_symbol == pdbstructure.extra_particle_identifier {
            ret.element = Some(EP.clone());
        } else {
            // try to find a sensible element symbol from columns 76-77
            let e = BY_SYMBOL.get(ret.element_symbol.as_str());
            if let Some(e) = e {
                ret.element = Some(e.clone());
            } else {
                ret.element = None;
            }
        }

        pdbstructure.next_atom_number = ret.serial_number + 1;
        pdbstructure.next_residue_number = ret.residue_number + 1;

        ret
    }
}

#[derive(Clone)]
struct Chain {
    chain_id: String,
    residues: Vec<()>,
    has_ter_record: bool,
    current_residue: Option<()>,
    residues_by_num_icode: HashMap<(), ()>,
    residues_by_number: HashMap<(), ()>,
}

impl Chain {
    fn new(chain_id: String) -> Self {
        Self {
            chain_id,
            residues: Vec::new(),
            has_ter_record: false,
            current_residue: None,
            residues_by_num_icode: HashMap::new(),
            residues_by_number: HashMap::new(),
        }
    }

    fn add_atom(&self, atom: Atom) {
        todo!()
    }
}

#[derive(Clone)]
struct Model {
    number: usize,
    chains: Vec<()>,
    current_chain: Option<Chain>,
    chains_by_id: HashMap<(), ()>,
    connect: Vec<()>,
}

impl Model {
    fn new(number: usize) -> Self {
        Self {
            number,
            chains: Vec::new(),
            current_chain: None,
            chains_by_id: HashMap::new(),
            connect: Vec::new(),
        }
    }

    fn add_atom(&mut self, atom: Atom) {
        if self.chains.is_empty() {
            self.add_chain(Chain::new(atom.chain_id.clone()));
        }
        self.add_chain(Chain::new(atom.chain_id.clone()));
        self.current_chain.as_mut().unwrap().add_atom(atom);
    }

    fn add_chain(&mut self, chain_id: Chain) {
        todo!()
    }
}

/// PdbStructure object holds a parsed Protein Data Bank format file.
///
/// Examples:
///
/// Load a pdb structure from a file:
/// > pdb = PdbStructure(open("1ARJ.pdb"))
///
/// Fetch the first atom of the structure:
/// > print pdb.iter_atoms().next()
/// ATOM      1  O5'   G N  17      13.768  -8.431  11.865  1.00  0.00           O
///
/// Loop over all of the atoms of the structure
/// > for atom in pdb.iter_atoms():
/// >     print atom
/// ATOM      1  O5'   G N  17      13.768  -8.431  11.865  1.00  0.00           O
/// ...
///
/// Get a list of all atoms in the structure:
/// > atoms = list(pdb.iter_atoms())
///
/// also:
/// residues = list(pdb.iter_residues())
/// positions = list(pdb.iter_positions())
/// chains = list(pdb.iter_chains())
/// models = list(pdb.iter_models())
///
/// Fetch atomic coordinates of first atom:
/// > print pdb.iter_positions().next()
/// [13.768, -8.431, 11.865] A
///
///  or
///
/// > print pdb.iter_atoms().next().position
/// [13.768, -8.431, 11.865] A
///
/// Strip the length units from an atomic position:
/// > import openmm.unit
/// > pos = pdb.iter_positions().next()
/// > print pos
/// [13.768, -8.431, 11.865] A
/// > print pos / openmm.unit.angstroms
/// [13.768, -8.431, 11.865]
/// > print pos / openmm.unit.nanometers
/// [1.3768, -0.8431, 1.1865]
///
///
/// The hierarchical structure of the parsed PDB structure is as follows:
/// PdbStructure
///   Model
///     Chain
///       Residue
///         Atom
///           Location
///
/// Model - A PDB structure consists of one or more Models.  Each model corresponds to one version of
/// an NMR structure, or to one frame of a molecular dynamics trajectory.
///
/// Chain - A Model contains one or more Chains.  Each chain corresponds to one molecule, although multiple
/// water molecules are frequently included in the same chain.
///
/// Residue - A Chain contains one or more Residues.  One Residue corresponds to one of the repeating
/// unit that constitutes a polymer such as protein or DNA.  For non-polymeric molecules, one Residue
/// represents one molecule.
///
/// Atom - A Residue contains one or more Atoms.  Atoms are chemical atoms.
///
/// Location - An atom can sometimes have more that one position, due to static disorder in X-ray
/// crystal structures.  To see all of the atom positions, use the atom.iter_positions() method,
/// or pass the parameter "include_alt_loc=True" to one of the other iter_positions() methods.
///
/// > for pos in pdb.iter_positions(include_alt_loc=True):
/// >   ...
///
/// Will loop over all atom positions, including multiple alternate locations for atoms that have
/// multiple positions.  The default value of include_alt_loc is False for the iter_positions()
/// methods.
#[derive(Default)]
struct PdbStructure {
    load_all_models: bool,
    extra_particle_identifier: String,
    models: Vec<Model>,
    current_model: Option<Model>,
    default_model: Option<Model>,
    models_by_number: HashMap<usize, Model>,
    periodic_box_vectors: Option<()>,
    sequences: Vec<()>,
    modified_residues: Vec<()>,
    next_atom_number: usize,
    next_residue_number: usize,
}

impl PdbStructure {
    pub fn new(
        filename: impl AsRef<Path>,
        load_all_models: bool,
        extra_particle_identifier: String,
    ) -> Result<Self, Box<dyn Error>> {
        let f = std::fs::File::open(filename)?;
        Ok(Self {
            load_all_models,
            extra_particle_identifier,
            models: todo!(),
            current_model: todo!(),
            default_model: todo!(),
            models_by_number: todo!(),
            periodic_box_vectors: todo!(),
            sequences: todo!(),
            modified_residues: todo!(),
            next_atom_number: 1,
            next_residue_number: 1,
        })
    }

    fn load(&mut self, f: impl Read) {
        let r = BufReader::new(f);
        let mut sel = Self::default();
        for pdb_line in r.lines().flatten() {
            let command = &pdb_line[..6];
            match command {
                "ATOM  " | "HETATM" => {
                    let atom = Atom::new(&pdb_line, self);
                    self.add_atom(atom);
                }
                _ => todo!("{command}"),
            }
        }
    }

    fn add_atom(&mut self, mut atom: Atom) {
        if self.current_model.is_none() {
            self.add_model(Model::new(0));
        }
        atom.model_number = self.current_model.as_ref().unwrap().number;
        self.current_model.as_mut().unwrap().add_atom(atom);
    }

    // TODO are these supposed to be references?? :O :(
    fn add_model(&mut self, model: Model) {
        if self.default_model.is_none() {
            self.default_model = Some(model.clone());
        }
        self.models.push(model.clone());
        self.current_model = Some(model.clone());
        if !self.models_by_number.contains_key(&model.number) {
            self.models_by_number.insert(model.number, model);
        }
    }
}

/// PDBFile parses a Protein Data Bank (PDB) file and constructs a Topology and
/// a set of atom positions from it.
///
/// This class also provides methods for creating PDB files. To write a file
/// containing a single model, call writeFile(). You also can create files that
/// contain multiple models. To do this, first call writeHeader(), then
/// writeModel() once for each model in the file, and finally writeFooter() to
/// complete the file.
pub struct PDBFile {
    topology: Topology,
    positions: Vec<Quantity<Vec3>>,
}

impl PDBFile {
    pub fn new(filename: impl AsRef<Path>) -> Self {
        let top = Topology::new();
        let pdb = PdbStructure::new(filename, true, "EP".to_owned());
        let positions = Vec::new();
        // TODO
        Self {
            topology: top,
            positions,
        }
    }
}

/// Modeller provides tools for editing molecular models, such as adding water
/// or missing hydrogens.
///
/// To use it, create a Modeller object, specifying the initial Topology and
/// atom positions. You can then call various methods to change the model in
/// different ways. Each time you do, a new Topology and list of coordinates is
/// created to represent the changed model. Finally, call getTopology() and
/// getPositions() to get the results.
pub struct Modeller {
    topology: Topology,
    positions: Vec<Quantity<Vec3>>,
}

impl Modeller {
    pub fn new(topology: Topology, positions: Vec<Quantity<Vec3>>) -> Self {
        Self {
            topology,
            positions,
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
pub struct Topology {
    chains: Vec<()>,
    num_residues: usize,
    num_atoms: usize,
    bonds: Vec<()>,
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
}

pub struct System;

pub struct Integrator;

/// Simulation provides a simplified API for running simulations with OpenMM and
/// reporting results.
///
/// A Simulation ties together various objects used for running a simulation: a
/// Topology, System, Integrator, and Context. To use it, you provide the
/// Topology, System, and Integrator, and it creates the Context automatically.
///
/// Simulation also maintains a list of "reporter" objects that record or
/// analyze data as the simulation runs, such as writing coordinates to files or
/// displaying structures on the screen. For example, the following line will
/// cause a file called "output.pdb" to be created, and a structure written to
/// it every 1000 time steps:
pub struct Simulation {
    topology: Topology,
    system: System,
    integrator: Integrator,
}

impl Simulation {
    pub fn new(
        topology: Topology,
        system: System,
        integrator: Integrator,
    ) -> Self {
        Self {
            topology,
            system,
            integrator,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// this is the part of the api I'm trying to replicate from ForceBalance
    #[test]
    fn openmmio() {
        let pdb = PDBFile::new("testfiles/test.pdb");
        let m = Modeller::new(pdb.topology, pdb.positions);
        let topology = m.topology;
        let system = todo!();
        let integrator = todo!();
        let simulation = Simulation::new(topology, system, integrator);
    }
}
