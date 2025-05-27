use derive_builder::Builder;

use crate::bindings::{self, CLevelSetExtractionParameters, CPoissonParams, CSolutionParameters};

/// This integer specifies the degree of the B-spline that is to be used to
/// define the finite elements system.
///
/// Larger degrees support higher order approximations, but come at the cost
/// of denser system matrices (incurring a cost in both space and time).
///
/// The default value for this parameter is 1.
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
#[repr(C)]
#[non_exhaustive]
pub enum Degree {
    /// Requires feature `degree_1`
    #[cfg(any(
        doc,
        feature = "degree_1",
        not(any(feature = "degree_1", feature = "degree_2", feature = "degree_3"))
    ))]
    One = 1,
    /// Requires feature `degree_2`
    #[cfg(any(doc, feature = "degree_2"))]
    Two = 2,
    /// Requires feature `degree_3`
    #[cfg(any(doc, feature = "degree_3"))]
    Three = 3,
}

/// This integer specifies the boundary type for the finite elements.
///
/// Valid values are:
/// 1: Free boundary constraints
/// 2: Dirichlet boundary constraints
/// 3: Neumann boundary constraints
///
/// The default value for this parameter is 3 (Neumann).
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
#[repr(C)]
#[non_exhaustive]
pub enum BoundaryType {
    /// Requires feature `boundary_free`
    #[cfg(any(doc, feature = "boundary_free"))]
    Free = 1,
    /// Requires feature `boundary_dirichlet`
    ///
    /// See <https://en.wikipedia.org/wiki/Dirichlet_boundary_condition>
    #[cfg(any(doc, feature = "boundary_dirichlet"))]
    Dirichlet = 2,
    /// Requires feature `boundary_neumann`
    ///
    /// See <https://en.wikipedia.org/wiki/Neumann_boundary_condition>
    #[cfg(any(
        doc,
        feature = "boundary_neumann",
        not(any(
            feature = "boundary_free",
            feature = "boundary_dirichlet",
            feature = "boundary_neumann"
        ))
    ))]
    Neumann = 3,
}

/// Options for Poisson Surface Reconstruction
#[derive(Builder, Debug, Clone, Default)]
#[builder(default, build_fn(private, name = "try_build"))]
pub struct PoissonParams {
    /// This integer specifies the degree of the B-spline that is to be used to
    /// define the finite elements system.
    ///
    /// Larger degrees support higher order approximations, but come at the cost
    /// of denser system matrices (incurring a cost in both space and time).
    ///
    /// The default value for this parameter is 1.
    pub degree: Degree,
    /// This integer specifies the boundary type for the finite elements.
    ///
    /// Valid values are:
    /// 1: Free boundary constraints
    /// 2: Dirichlet boundary constraints
    /// 3: Neumann boundary constraints
    ///
    /// The default value for this parameter is 3 (Neumann).
    pub boundary_type: BoundaryType,
    /// Main parameters for the Poisson surface reconstruction algorithm.
    pub solution_params: SolutionParameters,
    /// Parameters for the surface mesh extraction from the Poisson solution.
    pub level_set_extraction_params: LevelSetExtractionParameters,
}

/// Main parameters for the Poisson surface reconstruction algorithm.
#[derive(Builder, Debug, Clone, Default)]
#[builder(default, setter(strip_option), build_fn(private, name = "try_build"))]
pub struct SolutionParameters {
    verbose: Option<bool>,
    /// Defaults to false
    exact_interpolation: Option<bool>,
    /// Defaults to false
    show_residual: Option<bool>,
    /// Defaults to false
    confidence: Option<bool>,
    /// This floating point value specifies the ratio between the diameter of the
    /// cube used for reconstruction and the diameter of the samples' bounding cube.
    ///
    /// Defaults to 1.1
    scale: Option<f64>,
    /// Defaults to 0
    low_depth_cutoff: Option<f64>,
    /// This floating point value specifies the target width of the finest level octree cells.
    ///
    /// This parameter is ignored if the depth is also specified.
    width: Option<f64>,
    samples_per_node: Option<f64>,
    cg_solver_accuracy: Option<f64>,
    per_level_data_scale_factor: Option<f64>,
    /// This integer is the maximum depth of the tree that will be used for surface
    /// reconstruction. Running at depth d corresponds to solving on a grid whose
    /// resolution is no larger than 2^d x 2^d x ... Note that since the reconstructor
    /// adapts the octree to the sampling density, the specified reconstruction depth
    /// is only an upper bound.
    ///
    /// The default value for this parameter is 8.
    depth: Option<u32>,
    /// The depth up to which we can solve (solveDepth<=maxDepth).
    ///
    /// The finest level at which the system is solved.
    ///
    /// Defaults to depth.
    solve_depth: Option<u32>,
    /// The depth up to which the system is defined through the regular prolongation operators (baseDepth<=fullDepth).
    ///
    /// The coarsest level used in the multigrid solver.
    ///
    /// * A lower baseDepth relative to depth can sometimes produce smoother results.
    /// * This is because the finer details get interpolated from coarser approximations.
    ///
    /// Defaults to `full_depth`.
    base_depth: Option<u32>,
    /// The depth up to which the octree is completely refined (fullDepth<=maxDepth).
    ///
    /// Defaults to `solve_depth`.
    full_depth: Option<u32>,
    /// Controls density estimation resolution.
    ///
    /// # Practical effects:
    ///
    /// 1. For unevenly sampled points:
    /// * Higher kernelDepth helps the algorithm better adapt to varying point densities
    /// * Useful when some areas have dense sampling and others sparse sampling
    ///
    /// 2. For noisy data:
    /// * Lower kernelDepth can make the reconstruction more robust to noise
    /// * Creates a more averaged density estimate
    ///
    /// 3. For memory considerations:
    /// * Lower values use less memory during reconstruction
    /// * Can be important for very large point clouds
    ///
    /// 4. For reconstruction quality:
    /// * Setting to depth makes density estimation match the reconstruction resolution
    /// * Setting to 0 makes a uniform density estimate (less adaptive)
    ///
    /// Defaults to `depth-2` or 0.
    kernel_depth: Option<u32>,
    /// Defaults to 1
    base_vcycles: Option<u32>,
    /// This integer value specifies the number of Gauss-Seidel relaxations to be
    /// performed at each level of the hiearchy.
    ///
    /// The default value for this parameter is 8.
    iters: Option<u32>,
    /// Defaults to 0
    align_dir: Option<u32>,

    /// Defaults to false
    dirichlet_erode: Option<bool>,
    /// This floating point value specifies the importance that interpolation of
    /// the point samples is given in the formulation of the screened Poisson equation.
    /// The results of the original (unscreened) Poisson Reconstruction can be
    /// obtained by setting this value to 0. The default value for this parameter
    /// is twice the B-spline degree.
    point_weight: Option<f64>,
    /// Defaults to 0
    value_interpolation_weight: Option<f64>,
    envelope_depth: Option<u32>,
}

/// Parameters for the surface mesh extraction from the Poisson solution.
#[derive(Builder, Debug, Clone, Default)]
#[builder(default, setter(strip_option), build_fn(private, name = "try_build"))]
pub struct LevelSetExtractionParameters {
    /// Enabling this flag has the reconstructor use linear interpolation to estimate
    /// the positions of iso-vertices.
    ///
    /// Defaults to false.
    linear_fit: Option<bool>,
    /// Defaults to false.
    ///
    /// WARNING: NOT SUPPORTED BY THE BINDINGS YET.
    output_gradients: Option<bool>,
    /// Defaults to true.
    force_manifold: Option<bool>,
    /// Enabling this flag tells the reconstructor to output a polygon mesh (rather
    /// than triangulating the results of Marching Cubes).
    ///
    /// WARNING: NOT SUPPORTED BY THE BINDINGS YET.
    polygon_mesh: Option<bool>,
    /// Defaults to false.
    grid_coordinates: Option<bool>,
    /// Enabling this flag tells the reconstructor to output the estimated depth
    /// values of the iso-surface vertices.
    ///
    /// WARNING: NOT SUPPORTED BY THE BINDINGS YET.
    output_density: Option<bool>,
    verbose: Option<bool>,
}

impl SolutionParametersBuilder {
    pub fn build(&self) -> SolutionParameters {
        self.try_build().expect("to be infallible")
    }
}

impl SolutionParameters {
    fn apply(&self, params: &mut CSolutionParameters) {
        if let Some(v) = self.verbose {
            params._base.verbose = v;
        }
        if let Some(v) = self.exact_interpolation {
            params._base.exactInterpolation = v;
        }
        if let Some(v) = self.show_residual {
            params._base.showResidual = v;
        }
        if let Some(v) = self.confidence {
            params._base.confidence = v;
        }
        if let Some(v) = self.scale {
            params._base.scale = v;
        }
        if let Some(v) = self.low_depth_cutoff {
            params._base.lowDepthCutOff = v;
        }
        if let Some(v) = self.width {
            params._base.width = v;
        }
        if let Some(v) = self.samples_per_node {
            params._base.samplesPerNode = v;
        }
        if let Some(v) = self.cg_solver_accuracy {
            params._base.cgSolverAccuracy = v;
        }
        if let Some(v) = self.per_level_data_scale_factor {
            params._base.perLevelDataScaleFactor = v;
        }
        if let Some(v) = self.depth {
            params._base.depth = v;
        }
        if let Some(v) = self.solve_depth {
            params._base.solveDepth = v;
        }
        if let Some(v) = self.base_depth {
            params._base.baseDepth = v;
        }
        if let Some(v) = self.full_depth {
            params._base.fullDepth = v;
        }
        if let Some(v) = self.kernel_depth {
            params._base.kernelDepth = v;
        }
        if let Some(v) = self.base_vcycles {
            params._base.baseVCycles = v;
        }
        if let Some(v) = self.iters {
            params._base.iters = v;
        }
        if let Some(v) = self.align_dir {
            params._base.alignDir = v;
        }
        if let Some(v) = self.dirichlet_erode {
            params.dirichletErode = v;
        }
        if let Some(v) = self.point_weight {
            params.pointWeight = v;
        }
        if let Some(v) = self.value_interpolation_weight {
            params.valueInterpolationWeight = v;
        }
        if let Some(v) = self.envelope_depth {
            params.envelopeDepth = v;
        }
    }
}

impl LevelSetExtractionParametersBuilder {
    pub fn build(&self) -> LevelSetExtractionParameters {
        self.try_build().expect("to be infallible")
    }
}

impl PoissonParamsBuilder {
    pub fn build(&self) -> PoissonParams {
        self.try_build().expect("to be infallible")
    }
}

impl LevelSetExtractionParameters {
    fn apply(&self, params: &mut CLevelSetExtractionParameters) {
        if let Some(v) = self.linear_fit {
            params.linearFit = v;
        }
        if let Some(v) = self.output_gradients {
            params.outputGradients = v;
        }
        if let Some(v) = self.force_manifold {
            params.forceManifold = v;
        }
        if let Some(v) = self.polygon_mesh {
            params.polygonMesh = v;
        }
        if let Some(v) = self.grid_coordinates {
            params.gridCoordinates = v;
        }
        if let Some(v) = self.output_density {
            params.outputDensity = v;
        }
        if let Some(v) = self.verbose {
            params.verbose = v;
        }
    }
}

impl From<&PoissonParams> for CPoissonParams {
    fn from(override_params: &PoissonParams) -> Self {
        let mut params = unsafe { bindings::get_default_poisson_params() };

        override_params
            .solution_params
            .apply(&mut params.solution_params);
        override_params
            .level_set_extraction_params
            .apply(&mut params.level_set_extraction_params);

        params
    }
}

impl Default for Degree {
    fn default() -> Self {
        #[cfg(any(
            feature = "degree_1",
            not(any(feature = "degree_1", feature = "degree_2", feature = "degree_3"))
        ))]
        return Degree::One;

        #[cfg(all(not(feature = "degree_1"), feature = "degree_2"))]
        return Degree::Two;

        #[cfg(all(
            not(feature = "degree_1"),
            not(feature = "degree_2"),
            feature = "degree_3"
        ))]
        return Degree::Three;
    }
}

impl Default for BoundaryType {
    fn default() -> Self {
        #[cfg(any(
            feature = "boundary_neumann",
            not(any(
                feature = "boundary_free",
                feature = "boundary_dirichlet",
                feature = "boundary_neumann"
            ))
        ))]
        return BoundaryType::Neumann;

        #[cfg(all(not(feature = "boundary_neumann"), feature = "boundary_dirichlet"))]
        return BoundaryType::Dirichlet;

        #[cfg(all(
            not(feature = "boundary_neumann"),
            not(feature = "boundary_dirichlet"),
            feature = "boundary_free"
        ))]
        return BoundaryType::Free;
    }
}
