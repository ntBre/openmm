use std::ffi::c_int;

use openmm_sys::{
    OpenMM_State, OpenMM_State_destroy, OpenMM_State_getPositions,
    OpenMM_Vec3Array_get, OpenMM_Vec3Array_getSize,
};

use crate::{topology::Vec3, vec3};

bitflags::bitflags! {
    /// don't blame me, this is the name in C++
    #[repr(transparent)]
    pub struct DataType: u8 {
        const Positions = 1;
        const Velocities = 2;
        const Forces = 4;
        const Energy = 8;
        const Parameters = 16;
        const ParameterDerivatives = 32;
        const IntegratorParameters = 64;
    }
}

pub struct State {
    state: *mut OpenMM_State,
}

impl State {
    pub fn new(state: *mut OpenMM_State) -> Self {
        Self { state }
    }

    pub fn get_positions(&self) -> Vec<Vec3> {
        let positions =
            unsafe { OpenMM_State_getPositions(self.state as *const _) };
        let len = unsafe { OpenMM_Vec3Array_getSize(positions) };
        let mut ret = Vec::with_capacity(len as usize);
        for i in 0..len {
            let v = unsafe { OpenMM_Vec3Array_get(positions, i as c_int) };
            let (x, y, z) = unsafe { ((*v).x, (*v).y, (*v).z) };
            ret.push(vec3![x, y, z]);
        }
        ret
    }
}

impl Drop for State {
    fn drop(&mut self) {
        unsafe { OpenMM_State_destroy(self.state) }
    }
}
