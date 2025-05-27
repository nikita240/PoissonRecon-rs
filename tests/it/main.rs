//! This file is the entrypoint for the integration tests crate.
//!
//! Read [this blog](https://matklad.github.io/2021/02/27/delete-cargo-integration-tests.html)
//! to understand why we structure the tests this way.

use nalgebra::{Point3, Vector3};
use poissonrecon::{PoissonError, PoissonParams, PoissonParamsBuilder, SolutionParametersBuilder};

#[test]
fn test_empty_input() {
    let points = vec![];
    let normals = vec![];
    let result = poissonrecon::reconstruct_surface(
        &points,
        &normals,
        &PoissonParamsBuilder::default().build(),
    );

    assert!(result.is_err());
    assert_eq!(result.unwrap_err(), PoissonError::NoPoints);
}

#[test]
fn test_missing_normals() {
    let points = vec![Point3::new(0.0, 0.0, 0.0)];
    let normals = vec![];
    let result = poissonrecon::reconstruct_surface(
        &points,
        &normals,
        &PoissonParamsBuilder::default().build(),
    );

    assert!(result.is_err());
    assert_eq!(
        result.unwrap_err(),
        PoissonError::NormalsAndPointsMismatch(1, 0)
    );
}

#[test]
fn test_no_output_points() {
    let points = vec![Point3::new(0.0, 0.0, 0.0)];
    let normals = vec![Vector3::new(0.0, 0.0, 1.0)];
    let result = poissonrecon::reconstruct_surface(
        &points,
        &normals,
        &PoissonParamsBuilder::default().build(),
    );

    assert!(result.is_ok());
    let mesh = result.unwrap();
    assert!(mesh.vertices.is_empty());
    assert!(mesh.triangles.is_empty());
}

#[test]
fn test_min_valid_output() {
    let points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(-1.0, -1.0, 0.0),
    ];
    let normals = vec![
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.0, 0.0, 1.0),
    ];
    let result = poissonrecon::reconstruct_surface(
        &points,
        &normals,
        &PoissonParams {
            solution_params: SolutionParametersBuilder::default().depth(1).build(),
            ..Default::default()
        },
    );

    assert!(result.is_ok());
    let mesh = result.unwrap();
    assert!(!mesh.vertices.is_empty());
    assert!(!mesh.triangles.is_empty());
}
