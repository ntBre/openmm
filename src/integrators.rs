use std::ffi::c_int;

use openmm_sys::{
    OpenMM_VerletIntegrator, OpenMM_VerletIntegrator_create,
    OpenMM_VerletIntegrator_destroy, OpenMM_VerletIntegrator_step,
};

use crate::Integrator;

pub struct Verlet {
    integrator: *mut OpenMM_VerletIntegrator,
}

impl Verlet {
    pub fn new(step_size: f64) -> Self {
        let integrator = unsafe { OpenMM_VerletIntegrator_create(step_size) };
        Self { integrator }
    }

    pub fn step(&mut self, steps: i32) {
        unsafe {
            OpenMM_VerletIntegrator_step(self.integrator, steps as c_int);
        }
    }
}

impl Drop for Verlet {
    fn drop(&mut self) {
        unsafe {
            OpenMM_VerletIntegrator_destroy(self.integrator);
        }
    }
}

impl Integrator for Verlet {}
