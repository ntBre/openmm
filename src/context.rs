use bitflags::Flags;
use openmm_sys::{
    OpenMM_Context, OpenMM_Context_create, OpenMM_Context_destroy,
    OpenMM_Context_getState, OpenMM_Integrator, OpenMM_System,
};

use crate::{
    integrators::Integrator,
    state::{DataType, State},
    System,
};

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

    pub fn get_state(&self, types: DataType) -> State {
        // default value, might need to pass it
        let enforce_periodic_box = false;
        let state = unsafe {
            OpenMM_Context_getState(
                self.context as *const OpenMM_Context,
                types.bits() as i32,
                enforce_periodic_box as i32,
            )
        };
        State::new(state)
    }
}

impl Drop for Context {
    fn drop(&mut self) {
        unsafe {
            OpenMM_Context_destroy(self.context);
        }
    }
}
