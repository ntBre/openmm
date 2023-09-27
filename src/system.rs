use openmm_sys::{
    OpenMM_System, OpenMM_System_addConstraint, OpenMM_System_addParticle,
    OpenMM_System_create, OpenMM_System_destroy,
    OpenMM_System_getConstraintParameters, OpenMM_System_getNumConstraints,
    OpenMM_System_getNumForces, OpenMM_System_getNumParticles,
    OpenMM_System_getParticleMass, OpenMM_System_removeConstraint,
    OpenMM_System_removeForce, OpenMM_System_setConstraintParameters,
    OpenMM_System_setParticleMass,
    OpenMM_System_usesPeriodicBoundaryConditions,
};

pub struct System {
    pub(crate) system: *mut OpenMM_System,
}

impl System {
    pub fn new() -> Self {
        Self {
            system: unsafe { OpenMM_System_create() },
        }
    }

    /// get the number of particles in this System
    pub fn get_num_particles(&self) -> usize {
        unsafe { OpenMM_System_getNumParticles(self.system) as usize }
    }

    /// Add a particle to the System. If the mass is 0, Integrators will ignore
    /// the particle and not modify its position or velocity. This is most often
    /// used for virtual sites, but can also be used as a way to prevent a
    /// particle from moving. Returns the index of the particle that was added
    pub fn add_particle(&mut self, mass: f64) -> usize {
        unsafe { OpenMM_System_addParticle(self.system, mass) as usize }
    }

    pub fn get_particle_mass(&self, idx: usize) -> f64 {
        unsafe { OpenMM_System_getParticleMass(self.system, idx as i32) }
    }

    pub fn set_particle_mass(&mut self, idx: usize, mass: f64) {
        unsafe {
            OpenMM_System_setParticleMass(self.system, idx as i32, mass);
        }
    }

    // setVirtualSite
    // isVirtualSite
    // getVirtualSite

    /// get the number of distance constraints in this [System]
    pub fn get_num_constraints(&self) -> usize {
        unsafe { OpenMM_System_getNumConstraints(self.system) as usize }
    }

    /// Add a constraint to the System. Particles whose mass is 0 cannot
    /// participate in constraints. Returns the index of the constraint that was
    /// added
    pub fn add_constraint(
        &mut self,
        particle1: usize,
        particle2: usize,
        distance: f64,
    ) -> usize {
        unsafe {
            OpenMM_System_addConstraint(
                self.system,
                particle1 as i32,
                particle2 as i32,
                distance,
            ) as usize
        }
    }

    /// Get the parameters defining a distance constraint. Returns `particle1`
    /// (the index of the first particle involved in the constraint),
    /// `particle2`, and `distance` - the required distance between the two
    /// particles in nm
    pub fn get_constraint_parameters(&self, idx: usize) -> (usize, usize, f64) {
        let particle1 = std::ptr::null_mut();
        let particle2 = std::ptr::null_mut();
        let distance = std::ptr::null_mut();
        unsafe {
            OpenMM_System_getConstraintParameters(
                self.system,
                idx as i32,
                particle1,
                particle2,
                distance,
            );
        }
        (particle1 as usize, particle2 as usize, unsafe { *distance })
    }

    /// Set the parameters defining a distance constraint. Particles whose mass
    /// is 0 cannot participate in constraints.
    pub fn set_constraint_parameters(
        &mut self,
        idx: usize,
        particle1: usize,
        particle2: usize,
        distance: f64,
    ) {
        unsafe {
            OpenMM_System_setConstraintParameters(
                self.system,
                idx as i32,
                particle1 as i32,
                particle2 as i32,
                distance,
            );
        }
    }

    /// Remove a constraint from the [System].
    pub fn remove_constraint(&mut self, idx: usize) {
        unsafe {
            OpenMM_System_removeConstraint(self.system, idx as i32);
        }
    }

    // pub fn add_force(&mut self, force: Force) {}

    /// Get the number of Force objects that have been added to the System.
    pub fn get_num_forces(&self) -> usize {
        unsafe { OpenMM_System_getNumForces(self.system) as usize }
    }

    // pub fn get_force(&self, idx: usize) -> &Force {}

    // pub fn get_force_mut(&mut self, idx: usize) -> &mut Force {}

    /// Remove a Force from the System. The memory associated with the removed
    /// Force object is deleted.
    pub fn remove_force(&mut self, idx: usize) {
        unsafe {
            OpenMM_System_removeForce(self.system, idx as i32);
        }
    }

    // getDefaultPeriodicBoxVectors
    // setDefaultPeriodicBoxVectors

    /// Returns whether or not any forces in this System use periodic
    /// boundaries.
    pub fn uses_periodic_boundary_conditions(&self) -> bool {
        unsafe {
            OpenMM_System_usesPeriodicBoundaryConditions(self.system) != 0
        }
    }
}

impl Default for System {
    fn default() -> Self {
        Self::new()
    }
}

impl Drop for System {
    fn drop(&mut self) {
        unsafe {
            OpenMM_System_destroy(self.system);
        }
    }
}
