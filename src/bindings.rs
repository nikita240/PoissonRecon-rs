use std::ffi::c_int;

use crate::params::{BoundaryType, Degree};

type Real = std::ffi::c_double;

mod poisson_ctypes {
    #![allow(non_upper_case_globals)]
    #![allow(non_camel_case_types)]
    #![allow(non_snake_case)]
    #![allow(dead_code)]
    include!(concat!(env!("OUT_DIR"), "/poisson_ctypes.rs"));
}

pub type CSolutionParameters =
    poisson_ctypes::PoissonRecon_Reconstructor_Poisson_SolutionParameters<Real>;
pub type CLevelSetExtractionParameters =
    poisson_ctypes::PoissonRecon_Reconstructor_LevelSetExtractionParameters;

#[repr(C)]
pub struct CPoissonParams {
    pub solution_params:
        poisson_ctypes::PoissonRecon_Reconstructor_Poisson_SolutionParameters<Real>,
    pub level_set_extraction_params:
        poisson_ctypes::PoissonRecon_Reconstructor_LevelSetExtractionParameters,
}

unsafe extern "C" {
    pub fn poisson_reconstruct(
        data: *const PointNormalData,
        params: CPoissonParams,
        degree: Degree,
        boundary_type: BoundaryType,
    ) -> *mut MeshData;
    pub fn free_mesh_data(mesh: *mut MeshData);
    pub fn get_default_poisson_params() -> CPoissonParams;
}

#[repr(C)]
pub struct PointNormalData {
    pub points: *const [Real; 3],
    pub normals: *const [Real; 3],
    pub count: usize,
}

#[repr(C)]
pub struct MeshData {
    pub vertices: *mut Real,
    pub indices: *mut c_int,
    pub vertex_count: usize,
    pub index_count: usize,
}
