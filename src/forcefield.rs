use std::{collections::HashMap, error::Error, fs::read_to_string, path::Path};

use xmltree::Element;

use crate::{topology::Topology, System};

/// A ForceField constructs OpenMM System objects based on a Topology
#[derive(Default)]
pub struct ForceField {
    atom_types: HashMap<(), ()>,
    templates: HashMap<(), ()>,
    patches: HashMap<(), ()>,
    template_patches: HashMap<(), ()>,
    template_signatures: HashMap<(), ()>,
    atom_classes: HashMap<(), ()>,
    pub(crate) forces: Vec<()>,
    scripts: Vec<()>,
    template_matchers: Vec<()>,
    template_generators: Vec<()>,
}

impl ForceField {
    /// Load one (TODO: or more) XML files and create a ForceField object based
    /// on it (them).
    pub fn new(file: impl AsRef<Path>) -> Result<Self, Box<dyn Error>> {
        let mut ret = Self::default();
        ret.load_file(file)?;
        Ok(ret)
    }

    /// Load an XML file and add the definitions from it to this ForceField
    fn load_file(
        &mut self,
        file: impl AsRef<Path>,
    ) -> Result<(), Box<dyn Error>> {
        let data = read_to_string(file)?;
        let root = Element::parse(data.as_bytes())?;

        if let Some(types) = root.get_child("Include") {
            todo!();
        }

        if let Some(types) = root.get_child("AtomTypes") {
            todo!();
        }

        if let Some(types) = root.get_child("Residues") {
            todo!();
        }

        if let Some(types) = root.get_child("Patches") {
            todo!();
        }

        // TODO if child.tag in parsers. not sure how parsers gets initialized,
        // so far this doesn't seem to do anything at all

        if let Some(types) = root.get_child("Script") {
            todo!();
        }

        if let Some(types) = root.get_child("InitializationScript") {
            todo!();
        }

        Ok(())
    }

    /// Construct an OpenMM System representing a Topology with this force field
    pub fn create_system(&self, topology: Topology) -> System {
        todo!()
    }
}
