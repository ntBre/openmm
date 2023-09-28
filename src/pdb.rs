//! Load PDB files using the [PDBFile] struct. Internally, [PDBFile] creates a
//! [PdbStructure], which constructs a hierarchy of [Model] -> [Chain] ->
//! [Residue] -> [Atom]. By iterating along each of these, the [PDBFile] builds
//! its [Topology] and vector of atomic positions. A common pattern in this file
//! is the use of a `current_*` method that unwraps a `current_*` field on a
//! struct and indexes the field corresponding to the `*`. For example,
//! [PdbStructure::current_model] returns the
//! PdbStructure.current_model.unwrap()th element in `self.models`. This is in
//! contrast to the Python version, which simply stores a reference to the
//! current element in its self.current_* field

use crate::{
    element::{Element, BY_SYMBOL, EP},
    forcefield::ForceField,
    topology::{Topology, Vec3},
    Quantity,
};
use std::{
    collections::{HashMap, HashSet},
    error::Error,
    io::{BufRead, BufReader, Read},
    path::Path,
    sync::LazyLock,
};

#[allow(unused)]
#[derive(Clone, Debug)]
pub(crate) struct Location {
    pub(crate) alternate_location_indicator: String,
    pub(crate) position: Quantity<Vec3>,
    pub(crate) occupancy: f64,
    pub(crate) temperature_factor: Quantity<f64>,
    pub(crate) residue_name_with_spaces: String,
}

impl Location {
    pub(crate) fn new(
        alternate_location_indicator: String,
        xyz: Quantity<Vec3>,
        occupancy: f64,
        temperature_factor: Quantity<f64>,
        residue_name_with_spaces: String,
    ) -> Self {
        Self {
            alternate_location_indicator,
            position: xyz,
            occupancy,
            temperature_factor,
            residue_name_with_spaces,
        }
    }
}

/// Atom represents one atom in a PDB structure.
#[derive(Clone, Debug, Default)]
pub(crate) struct Atom {
    pub(crate) is_first_atom_in_chain: bool,
    pub(crate) is_final_atom_in_chain: bool,
    pub(crate) is_first_residue_in_chain: bool,
    pub(crate) is_final_residue_in_chain: bool,
    #[allow(unused)]
    pub(crate) record_name: String,
    pub(crate) serial_number: usize,
    pub(crate) name_with_spaces: String,
    pub(crate) residue_name_with_spaces: String,
    pub(crate) residue_name: String,
    pub(crate) chain_id: String,
    pub(crate) residue_number: usize,
    pub(crate) insertion_code: String,
    pub(crate) locations: HashMap<String, Location>,
    pub(crate) default_location_id: String,
    pub(crate) segment_id: String,
    pub(crate) element_symbol: String,
    pub(crate) formal_charge: Option<isize>,
    pub(crate) element: Option<Element>,
    pub(crate) model_number: usize,
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
    pub(crate) fn new(pdb_line: &str, pdbstructure: &mut PdbStructure) -> Self {
        let mut ret = Self {
            record_name: pdb_line[..6].trim().to_owned(),
            // TODO this could also be in hex supposedly
            serial_number: match pdb_line[6..11].parse() {
                Ok(i) => i,
                Err(_) => pdbstructure.next_atom_number,
            },
            name_with_spaces: pdb_line[12..16].to_owned(),
            ..Self::default()
        };

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

        ret.residue_number = match pdb_line[22..26].trim().parse() {
            Ok(n) => n,
            Err(e) => panic!(
                "failed to parse {} as res number with {e}",
                &pdb_line[22..26]
            ),
        };

        ret.insertion_code = pdb_line[26..27].to_owned();

        let x: f64 = pdb_line[30..38]
            .trim()
            .parse()
            .expect("failed to parse x coord");
        let y: f64 = pdb_line[38..46]
            .trim()
            .parse()
            .expect("failed to parse y coord");
        let z: f64 = pdb_line[46..54]
            .trim()
            .parse()
            .expect("failed to parse z coord");

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

    pub(crate) fn alternate_location_indicator(&self) -> &str {
        let id = &self.default_location_id;
        self.locations[id.as_str()]
            .alternate_location_indicator
            .as_str()
    }

    pub(crate) fn get_name(&self) -> &String {
        &self.residue_name
    }

    pub(crate) fn get_position(&self) -> &Quantity<Vec3> {
        &self.get_location().position
    }

    pub(crate) fn get_location(&self) -> &Location {
        let id = &self.default_location_id;
        &self.locations[id]
    }
}

// TODO could just be Location if Resiude has its own module. internal to
// Residue definition in Python
#[allow(unused)]
#[derive(Clone, Debug)]
pub(crate) struct ResLoc {
    pub(crate) alternate_location_indicator: String,
    pub(crate) residue_name_with_spaces: String,
}

impl ResLoc {
    pub(crate) fn new(
        alternate_location_indicator: String,
        residue_name_with_spaces: String,
    ) -> Self {
        Self {
            alternate_location_indicator,
            residue_name_with_spaces,
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Residue {
    #[allow(unused)]
    pub(crate) primary_location_id: String,
    pub(crate) locations: HashMap<String, ResLoc>,
    pub(crate) name_with_spaces: String,
    pub(crate) number: usize,
    pub(crate) insertion_code: String,
    pub(crate) atoms: Vec<Atom>,
    pub(crate) atoms_by_name: HashMap<String, Atom>,
    pub(crate) is_first_in_chain: bool,
    pub(crate) is_final_in_chain: bool,
    pub(crate) current_atom: Option<usize>,
}

impl Residue {
    pub(crate) fn new(
        name: String,
        number: usize,
        insertion_code: String, // default = ' '
        primary_alternate_location_indicator: String, // default = ' '
    ) -> Self {
        let alt_loc = primary_alternate_location_indicator;
        Self {
            primary_location_id: alt_loc.clone(),
            // wtf is this?? is it supposed to be a reference? why do you store
            // alt_loc in the value of the map if the key is also the alt_loc??
            locations: HashMap::from([(
                alt_loc.clone(),
                ResLoc::new(alt_loc, name.clone()),
            )]),
            name_with_spaces: name,
            number,
            insertion_code,
            atoms: Vec::new(),
            atoms_by_name: HashMap::new(),
            is_first_in_chain: false,
            is_final_in_chain: false,
            current_atom: None,
        }
    }

    pub(crate) fn get_name(&self) -> &str {
        let alt_loc = &self.primary_location_id;
        let loc = &self.locations[alt_loc];
        loc.residue_name_with_spaces.trim()
    }

    pub(crate) fn finalize(&mut self) {
        if !self.atoms.is_empty() {
            self.atoms[0].is_first_atom_in_chain = self.is_first_in_chain;
            self.atoms.last_mut().unwrap().is_final_atom_in_chain =
                self.is_final_in_chain;
            for atom in self.atoms.iter_mut() {
                atom.is_first_residue_in_chain = self.is_first_in_chain;
                atom.is_final_residue_in_chain = self.is_final_in_chain;
            }
        }
    }

    pub(crate) fn add_atom(&mut self, atom: Atom) {
        let alt_loc = atom.alternate_location_indicator();
        if !self.locations.contains_key(alt_loc) {
            self.locations.insert(
                alt_loc.to_owned(),
                ResLoc::new(
                    alt_loc.to_owned(),
                    atom.residue_name_with_spaces.clone(),
                ),
            );
        }
        assert_eq!(atom.residue_number, self.number);
        assert_eq!(atom.insertion_code, self.insertion_code);
        if self.atoms_by_name.contains_key(&atom.name_with_spaces) {
            let old_atom =
                self.atoms_by_name.get_mut(&atom.name_with_spaces).unwrap();
            if old_atom
                .locations
                .contains_key(atom.alternate_location_indicator())
            {
                eprintln!("warning duplicate atoms");
            } else {
                for (alt_loc, position) in atom.locations {
                    old_atom.locations.insert(alt_loc, position);
                }
                return; // no new atom added
            }
        }
        self.atoms_by_name
            .insert(atom.get_name().to_owned(), atom.clone());
        self.atoms_by_name
            .insert(atom.name_with_spaces.clone(), atom.clone());
        self.atoms.push(atom);
        self.current_atom = Some(self.atoms.len() - 1);
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Chain {
    pub(crate) chain_id: String,
    pub(crate) residues: Vec<Residue>,
    pub(crate) has_ter_record: bool,
    pub(crate) current_residue: Option<usize>,
    pub(crate) residues_by_num_icode: HashMap<String, Residue>,
    pub(crate) residues_by_number: HashMap<usize, Residue>,
}

impl Chain {
    pub(crate) fn new(chain_id: String) -> Self {
        Self {
            chain_id,
            residues: Vec::new(),
            has_ter_record: false,
            current_residue: None,
            residues_by_num_icode: HashMap::new(),
            residues_by_number: HashMap::new(),
        }
    }

    pub(crate) fn add_atom(&mut self, atom: Atom) {
        if self.residues.is_empty() {
            let residue = Residue::new(
                atom.residue_name_with_spaces.clone(),
                atom.residue_number,
                atom.insertion_code.clone(),
                atom.alternate_location_indicator().to_owned(),
            );
            self.add_residue(residue);
        } else if self.current_residue().number != atom.residue_number {
            panic!("case 2");
        } else if self.current_residue().insertion_code != atom.insertion_code {
            panic!("case 3");
        } else if self.current_residue().name_with_spaces
            != atom.residue_name_with_spaces
        {
            panic!("case 4");
        } else if atom.alternate_location_indicator() != " " {
            panic!("case 5");
        } else {
            //warning
            eprintln!("warning: two consecutive residues with the same number");
        }
        self.current_residue().add_atom(atom);
    }

    pub(crate) fn current_residue(&mut self) -> &mut Residue {
        &mut self.residues[self.current_residue.unwrap()]
    }

    pub(crate) fn iter_residues(&self) -> std::slice::Iter<'_, Residue> {
        self.residues.iter()
    }

    pub(crate) fn add_residue(&mut self, mut residue: Residue) {
        if self.residues.is_empty() {
            residue.is_first_in_chain = true;
        }
        self.residues.push(residue.clone());
        self.current_residue = Some(self.residues.len() - 1);
        let key = residue.number.to_string() + &residue.insertion_code;
        self.residues_by_num_icode
            .entry(key)
            .or_insert(residue.clone());
        self.residues_by_number
            .entry(residue.number)
            .or_insert(residue);
    }

    pub(crate) fn add_ter_record(&mut self) {
        self.has_ter_record = true;
        self.finalize();
    }

    pub(crate) fn finalize(&mut self) {
        self.residues[0].is_first_in_chain = true;
        self.residues.last_mut().unwrap().is_final_in_chain = true;
        for residue in self.residues.iter_mut() {
            residue.finalize();
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Model {
    pub(crate) number: usize,
    pub(crate) chains: Vec<Chain>,
    pub(crate) current_chain: Option<usize>,
    pub(crate) chains_by_id: HashMap<String, Chain>,
    pub(crate) connects: Vec<Vec<usize>>,
}

impl Model {
    pub(crate) fn new(number: usize) -> Self {
        Self {
            number,
            chains: Vec::new(),
            current_chain: None,
            chains_by_id: HashMap::new(),
            connects: Vec::new(),
        }
    }

    pub(crate) fn current_chain(&mut self) -> &mut Chain {
        &mut self.chains[self.current_chain.unwrap()]
    }

    pub(crate) fn add_atom(&mut self, atom: Atom) {
        if self.chains.is_empty() {
            self.add_chain(Chain::new(atom.chain_id.clone()));
        }
        // shouldn't these be else if with above?? or can only one be true at
        // once? either way else if would clarify the intent
        if self.current_chain().chain_id != atom.chain_id
            || self.current_chain().has_ter_record
        {
            self.add_chain(Chain::new(atom.chain_id.clone()));
        }
        self.current_chain().add_atom(atom);
    }

    pub(crate) fn add_chain(&mut self, chain: Chain) {
        self.chains.push(chain.clone());
        self.current_chain = Some(self.chains.len() - 1);
        self.chains_by_id
            .entry(chain.chain_id.clone())
            .or_insert(chain);
    }

    pub(crate) fn finalize(&mut self) {
        for chain in self.chains.iter_mut() {
            chain.finalize();
        }
    }

    pub(crate) fn iter_chains(&self) -> std::slice::Iter<'_, Chain> {
        self.chains.iter()
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
#[allow(unused)]
#[derive(Default)]
pub(crate) struct PdbStructure {
    pub(crate) load_all_models: bool,
    pub(crate) extra_particle_identifier: String,
    pub(crate) models: Vec<Model>,
    pub(crate) current_model: Option<usize>,
    pub(crate) default_model: Option<Model>,
    pub(crate) models_by_number: HashMap<usize, Model>,
    pub(crate) periodic_box_vectors: Option<()>,
    pub(crate) sequences: Vec<()>,
    pub(crate) modified_residues: Vec<()>,
    pub(crate) next_atom_number: usize,
    pub(crate) next_residue_number: usize,
}

impl PdbStructure {
    pub fn new(
        filename: impl AsRef<Path>,
        load_all_models: bool,
        extra_particle_identifier: String,
    ) -> Result<Self, Box<dyn Error>> {
        let f = std::fs::File::open(filename)?;
        let mut ret = Self {
            load_all_models,
            extra_particle_identifier,
            models: Vec::new(),
            current_model: None,
            default_model: None,
            models_by_number: HashMap::new(),
            periodic_box_vectors: None,
            sequences: Vec::new(),
            modified_residues: Vec::new(),
            next_atom_number: 1,
            next_residue_number: 1,
        };
        ret.load(f);
        Ok(ret)
    }

    pub(crate) fn current_model(&mut self) -> &mut Model {
        &mut self.models[self.current_model.unwrap()]
    }

    pub(crate) fn load(&mut self, f: impl Read) {
        self.reset_atom_numbers();
        self.reset_residue_numbers();
        let r = BufReader::new(f);
        for pdb_line in r.lines().flatten() {
            let command = &pdb_line[..6];
            match command {
                "ATOM  " | "HETATM" => {
                    let atom = Atom::new(&pdb_line, self);
                    self.add_atom(atom);
                }
                "CONECT" => {
                    let atoms = pdb_line[6..]
                        .split_ascii_whitespace()
                        .flat_map(|s| s.parse::<usize>())
                        .collect();
                    self.current_model().connects.push(atoms);
                }
                "MODEL " => todo!(),
                "ENDMDL" => todo!(),
                "END   " => todo!(),
                "TER   " => {
                    self.current_model().current_chain().add_ter_record();
                    self.reset_residue_numbers();
                }
                "CRYST1" => todo!(),
                "SEQRES" => todo!(),
                "MODRES" => todo!(),
                _ => {}
            }
        }
        self.finalize();
    }

    pub(crate) fn add_atom(&mut self, mut atom: Atom) {
        if self.current_model.is_none() {
            self.add_model(Model::new(0));
        }
        atom.model_number = self.current_model().number;
        self.current_model().add_atom(atom);
    }

    // TODO are these supposed to be references?? :O :(
    pub(crate) fn add_model(&mut self, model: Model) {
        if self.default_model.is_none() {
            self.default_model = Some(model.clone());
        }
        self.models.push(model.clone());
        self.current_model = Some(self.models.len() - 1);
        self.models_by_number.entry(model.number).or_insert(model);
    }

    // TODO this should return an iterator, but I can't figure out how to get
    // the signature right
    pub(crate) fn iter_chains(&self) -> Vec<Chain> {
        self.models
            .iter()
            .flat_map(|model| model.chains.clone().into_iter())
            .collect()
    }

    pub(crate) fn reset_atom_numbers(&mut self) {
        self.next_atom_number = 1;
    }

    pub(crate) fn reset_residue_numbers(&mut self) {
        self.next_residue_number = 1;
    }

    pub(crate) fn finalize(&mut self) {
        for model in self.models.iter_mut() {
            model.finalize();
        }
    }

    /// the default for `use_all_models` is actually false, in which case only
    /// the first model is returned
    pub(crate) fn iter_models(
        &self,
        use_all_models: bool,
    ) -> Box<dyn Iterator<Item = &Model> + '_> {
        if use_all_models {
            Box::new(self.models.iter())
        } else {
            Box::new(self.models.iter().take(1))
        }
    }
}

// TODO actually build this thing
pub(crate) static RESIDUE_NAME_REPLACEMENTS: LazyLock<HashMap<String, String>> =
    LazyLock::new(HashMap::new);

/// PDBFile parses a Protein Data Bank (PDB) file and constructs a Topology and
/// a set of atom positions from it.
///
/// This class also provides methods for creating PDB files. To write a file
/// containing a single model, call writeFile(). You also can create files that
/// contain multiple models. To do this, first call writeHeader(), then
/// writeModel() once for each model in the file, and finally writeFooter() to
/// complete the file.
#[derive(Debug)]
pub struct PDBFile {
    pub topology: Topology,
    pub positions: Vec<Quantity<Vec<Vec3>>>,
}

impl PDBFile {
    pub fn new(filename: impl AsRef<Path>) -> Self {
        let mut top = Topology::new();
        let pdb = PdbStructure::new(filename, true, "EP".to_owned()).unwrap();

        // build the topology
        for chain in pdb.iter_chains() {
            let c = top.add_chain(chain.chain_id.clone());
            for residue in chain.iter_residues() {
                let mut res_name = residue.get_name();
                if let Some(name) = RESIDUE_NAME_REPLACEMENTS.get(res_name) {
                    res_name = name;
                }
                let _r = top.add_residue(
                    res_name,
                    c.clone(),
                    residue.number.to_string(),
                    residue.insertion_code.clone(),
                );
            }
        }
        let mut positions = Vec::new();
        for model in pdb.iter_models(true) {
            let mut coords = Vec::new();
            for chain in model.iter_chains() {
                for residue in chain.iter_residues() {
                    let mut processed_atom_names = HashSet::new();
                    for atom in residue.atoms_by_name.values() {
                        if processed_atom_names.contains(atom.get_name())
                            || atom.residue_name != residue.get_name()
                        {
                            continue;
                        }
                        processed_atom_names.insert(atom.get_name());
                        // TODO ensure nm
                        let Quantity { value: pos, .. } = atom.get_position();
                        coords.push(pos.clone());
                    }
                }
            }
            positions.push(Quantity::new(coords, "nanometers"));
        }

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
    pub topology: Topology,
    pub positions: Vec<Quantity<Vec<Vec3>>>,
}

impl Modeller {
    pub fn new(
        topology: Topology,
        positions: Vec<Quantity<Vec<Vec3>>>,
    ) -> Self {
        Self {
            topology,
            positions,
        }
    }

    /// Add missing extra particles to the model that are required by the force
    /// field.
    ///
    /// Some force fields use "extra particles" that do not represent actual
    /// atoms, but still need to be included in the System. Examples include
    /// lone pairs, Drude particles, and the virtual sites used in some water
    /// models to adjust the charge distribution. Extra particles can be
    /// recognized by the fact that their element is None.
    ///
    /// This method is primarily used to add extra particles, but it can also
    /// remove them (...). It tries to match every residue in the Topology to a
    /// template in the force field. If there is no match, it will both add and
    /// remove extra particles as necessary to make it match.
    #[allow(unused)]
    pub(crate) fn add_extra_particles(&mut self, forcefield: ForceField) {
        // record which atoms are bonded to each other atom
        let mut bonded_to_atom =
            vec![HashSet::new(); self.topology.atoms().count()];
        for (atom1, atom2) in self.topology.bonds() {
            bonded_to_atom[atom1.index].insert(atom2.index);
            bonded_to_atom[atom2.index].insert(atom1.index);
        }

        // if the force field has a DrudeForce, record the types of Drude
        // particles and their parents since we'll need them for picking
        // positions
        for _force in forcefield.forces {
            todo!();
        }

        // identify the template to use for each residue

        todo!()
    }

    pub fn get_topology(&self) -> &Topology {
        &self.topology
    }
}
