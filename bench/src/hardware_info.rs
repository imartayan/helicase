use sysinfo::System;

#[cfg(target_feature = "avx2")]
pub const VECTOR_TECH: &str = "avx2";
#[cfg(target_feature = "neon")]
pub const VECTOR_TECH: &str = "neon";
#[cfg(not(any(target_feature = "avx2", target_feature = "neon")))]
pub const VECTOR_TECH: &str = "none";

pub struct CpuInfo {
    pub brand: String,
    pub vendor_id: String,
    pub vector_tech: String,
}

pub fn get_hardware_info() -> CpuInfo {
    let mut sys = System::new();
    sys.refresh_cpu_all();

    let cpu = &sys.cpus()[0];
    CpuInfo {
        brand: cpu.brand().to_string(),
        vendor_id: cpu.vendor_id().to_string(),
        vector_tech: VECTOR_TECH.to_string(),
    }
}
