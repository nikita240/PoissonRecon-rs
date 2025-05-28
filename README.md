[![CI](https://github.com/nikita240/PoissonRecon-rs/actions/workflows/ci.yaml/badge.svg)](https://github.com/nikita240/PoissonRecon-rs/actions/workflows/ci.yaml) [![crates.io](https://img.shields.io/crates/v/poissonrecon.svg)](https://crates.io/crates/poissonrecon) [![docs.rs](https://img.shields.io/docsrs/poissonrecon)](https://docs.rs/poissonrecon)

# PoissonRecon-rs

<!-- cargo-rdme start -->

Safe Rust bindings for the Poisson Surface Reconstruction algorithm implemented
by Michael Kazhdan.

Original source: <https://github.com/mkazhdan/PoissonRecon>

## Example

```rust
let mesh = poissonrecon::reconstruct_surface(
    &points,
    &normals,
    &poissonrecon::PoissonParamsBuilder::default().build(),
)?;
```

## Comparison to the [`poisson_reconstruction`](https://crates.io/crates/poisson_reconstruction) crate

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

## Template unrolling

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

<!-- cargo-rdme end -->
