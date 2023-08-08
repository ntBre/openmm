use openmm_sys::{OpenMM_System, OpenMM_System_create, OpenMM_System_destroy};

pub struct System {
    pub(crate) system: *mut OpenMM_System,
}

impl System {
    pub fn new() -> Self {
        Self {
            system: unsafe { OpenMM_System_create() },
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
