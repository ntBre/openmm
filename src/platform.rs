use std::ffi::CString;

use openmm_sys::{
    OpenMM_Platform, OpenMM_Platform_destroy, OpenMM_Platform_getPlatformByName,
};

pub struct Platform {
    platform: *mut OpenMM_Platform,
}

impl Platform {
    pub fn by_name(name: &str) -> Self {
        let name = CString::new(name).unwrap();
        let platform =
            unsafe { OpenMM_Platform_getPlatformByName(name.as_ptr()) };
        Self { platform }
    }
}

impl Drop for Platform {
    fn drop(&mut self) {
        unsafe { OpenMM_Platform_destroy(self.platform) }
    }
}
