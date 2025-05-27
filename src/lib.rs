/*!
Safe Rust bindings for the Poisson Surface Reconstruction algorithm implemented
by Michael Kazhdan.

Original source: <https://github.com/mkazhdan/PoissonRecon>

# Example

```
# fn main() -> Result<(), poissonrecon::PoissonError> {
#     let points = vec![nalgebra::Point3::new(0.0, 0.0, 0.0)];
#     let normals = vec![nalgebra::Vector3::new(0.0, 0.0, 1.0)];
let mesh = poissonrecon::reconstruct_surface(
    &points,
    &normals,
    &poissonrecon::PoissonParamsBuilder::default().build(),
)?;
#     Ok(())
# }
```

# Comparison to the [`poisson_reconstruction`](https://crates.io/crates/poisson_reconstruction) crate

The `poisson_reconstruction` crate is an excellent Rust reimplementation of
the Screened Poisson Surface Reconstruction algorithm. It's maintained by
the Foresight Mining Software Corporation.

This `poissonrecon` crate on the other hand, is a Rust wrapper around Michael Kazhdan's
original C++ implementation. It compiles a C++ library that you must drag along
with your Rust project. It also exposes more tuning parameters for the algorithm
than the `poisson_reconstruction` crate does, which may be useful for some
applications.

You should use the `poisson_reconstruction` crate if you want a pure Rust
implementation.

You should use the `poissonrecon` crate if you want to use the original battle-tested
implementation, and you don't mind dragging along a C++ library.

The repository contains an example comparing the output of the two implementations:

```sh
cargo run --release --example comparison
```

# Template unrolling

The PoissonRecon C++ algorithm is templated over the FEMDegree and the BoundaryType.
You can enable compilation of different combinations of FEMDegree and BoundaryType
by enabling or disabling the corresponding feature flags.

These are the available feature flags:
- `degree_1`: Enables the [`Degree::One`] option.
- `degree_2`: Enables the [`Degree::Two`] option.
- `degree_3`: Enables the [`Degree::Three`] option.
- `boundary_free`: Enables the [`BoundaryType::Free`] option.
- `boundary_dirichlet`: Enables the [`BoundaryType::Dirichlet`] option.
- `boundary_neumann`: Enables the [`BoundaryType::Neumann`] option.

`degree_1` and `boundary_neumann` are enabled by default.

> ⚠️ You must have at least one degree and one boundary type enabled. If none
> are enabled, the library will forcefully compile as if `degree_1` and `boundary_neumann`
> were enabled.

It is recommended to enable only the feature flags you need, as each
combination will significantly increase the compilation time and binary size.
Enabling all feature flags will result in a `3 * 3 = 9` increase in compilation time
and binary size.
*/

use std::slice;

use bytemuck::cast_slice;
use nalgebra::{Point3, Vector3};
use thiserror::Error;

pub(crate) mod bindings;
mod params;

pub use crate::params::{
    BoundaryType, Degree, LevelSetExtractionParameters, LevelSetExtractionParametersBuilder,
    PoissonParams, PoissonParamsBuilder, SolutionParameters, SolutionParametersBuilder,
};

/// Perform Poisson Surface Reconstruction on the given points and normals.
///
/// # Arguments
///
/// * `points` - The input points to reconstruct from.
/// * `normals` - The normals corresponding to each point.
/// * `params` - Poisson configuration.
///
/// # Returns
///
/// A mesh containing vertices and triangle indices.
pub fn reconstruct_surface(
    points: &[Point3<f64>],
    normals: &[Vector3<f64>],
    params: &PoissonParams,
) -> Result<Mesh3D, PoissonError> {
    if points.is_empty() {
        return Err(PoissonError::NoPoints);
    }

    if points.len() != normals.len() {
        return Err(PoissonError::NormalsAndPointsMismatch(
            points.len(),
            normals.len(),
        ));
    }

    // Prepare input data
    // Make sure this doesn't get dropped until poisson is done!!!
    let data = bindings::PointNormalData {
        points: cast_slice(points).as_ptr(),
        normals: cast_slice(normals).as_ptr(),
        count: points.len(),
    };

    // # Safety
    //
    // This whole block is unsafe because we are letting C++ code operate on data
    // allocated by Rust, after which we operate on data allocated by C++ and
    // then manually call a C function to free the C++ allocated memory.
    unsafe {
        // # Safety
        //
        // MeshData is allocated entirely by C++ so we cannot manage it's memory
        // with Rust. We intentionally leak it from RAII C++ objects so that
        // we can let Rust copy out the data and then free it manually at the end
        // with `bindings::free_mesh_data()`.
        let mesh_data_ptr = bindings::poisson_reconstruct(
            &data,
            params.into(),
            params.degree,
            params.boundary_type,
        );
        if mesh_data_ptr.is_null() {
            return Err(PoissonError::ReconstructionFailed);
        }

        let mesh_data = &*mesh_data_ptr;

        let mut vertices = Vec::with_capacity(mesh_data.vertex_count);
        let mut indices = Vec::with_capacity(mesh_data.index_count);

        // Convert C mesh to Rust Mesh
        //
        // # Safety
        //
        // The `.vertices` and `.indices` may be null if no data is returned. We
        // check for that.
        //
        // The `mesh_data.vertex_count` and `mesh_data.index_count` represent the
        // number of vertices and triangles, while the raw data itself is a flat
        // array of numbers, which is why we multiply the count by 3.
        {
            if !mesh_data.vertices.is_null() && mesh_data.vertex_count != 0 {
                let vertices_slice =
                    slice::from_raw_parts(mesh_data.vertices, mesh_data.vertex_count * 3);

                // # Safety
                //
                // We HAVE to eat one copy because the structures on the C side are
                // not allocated by the same allocator, so we can't take ownership of them.
                //
                // See [`Vec::from_raw_parts`](https://doc.rust-lang.org/std/vec/struct.Vec.html#method.from_raw_parts)
                // for more details on why we can't use the data directly.
                for i in 0..mesh_data.vertex_count {
                    vertices.push(Point3::new(
                        vertices_slice[i * 3],
                        vertices_slice[i * 3 + 1],
                        vertices_slice[i * 3 + 2],
                    ));
                }
            }
            if !mesh_data.indices.is_null() && mesh_data.index_count != 0 {
                let indices_slice =
                    slice::from_raw_parts(mesh_data.indices, mesh_data.index_count * 3);

                // Resize it at the same time as we're copying it.
                //
                // # Safety
                //
                // We HAVE to eat one copy because the structures on the C side are
                // not allocated by the same allocator, so we can't take ownership of them.
                //
                // See [`Vec::from_raw_parts`](https://doc.rust-lang.org/std/vec/struct.Vec.html#method.from_raw_parts)
                // for more details on why we can't use the data directly.
                for i in 0..mesh_data.index_count {
                    indices.push([
                        indices_slice[i * 3] as usize,
                        indices_slice[i * 3 + 1] as usize,
                        indices_slice[i * 3 + 2] as usize,
                    ]);
                }
            }
        }

        // Free C allocated memory
        bindings::free_mesh_data(mesh_data_ptr);

        Ok(Mesh3D {
            vertices,
            triangles: indices,
        })
    }
}

/// Simple 3D mesh datastructure.
#[derive(Clone, Default, Debug)]
pub struct Mesh3D {
    /// Position of each vertex.
    pub vertices: Vec<Point3<f64>>,

    /// The mesh triangles.
    ///
    /// Each entry in this array contains the indexes of the 3 vertices that make
    /// up the triangle.
    pub triangles: Vec<[usize; 3]>,
}

#[derive(Error, Debug, PartialEq, Eq)]
pub enum PoissonError {
    #[error("No points provided")]
    NoPoints,
    #[error("Got {0} points and {1} normals, but they must match")]
    NormalsAndPointsMismatch(usize, usize),
    #[error("PoissonRecon C++ solver returned null")]
    ReconstructionFailed,
}
