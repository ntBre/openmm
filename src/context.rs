use openmm_sys::{
    OpenMM_Context, OpenMM_Context_create, OpenMM_Context_destroy,
    OpenMM_Integrator, OpenMM_System,
};

use crate::{integrators::Integrator, System};

pub struct Context {
    context: *mut OpenMM_Context,
}

impl Context {
    pub fn new<I: Integrator>(system: &System, integrator: &mut I) -> Self {
        Self {
            context: unsafe {
                OpenMM_Context_create(system.system, integrator.integrator())
            },
        }
    }
}

impl Drop for Context {
    fn drop(&mut self) {
        unsafe {
            OpenMM_Context_destroy(self.context);
        }
    }
}
