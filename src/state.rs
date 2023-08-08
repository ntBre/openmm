use openmm_sys::OpenMM_State;

/// don't blame me, this is the name in C++
bitflags::bitflags! {
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
}
