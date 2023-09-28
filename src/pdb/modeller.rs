use std::collections::HashSet;

use crate::forcefield::ForceField;
use crate::quantity::Quantity;
use crate::topology::Topology;
use crate::topology::Vec3;

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
    pub fn add_extra_particles(&mut self, forcefield: ForceField) {
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
