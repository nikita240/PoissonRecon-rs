#![allow(deprecated)]

use bevy::pbr::wireframe::WireframePlugin;
use bevy::prelude::*;
use bevy_panorbit_camera::{PanOrbitCamera, PanOrbitCameraPlugin};
use nalgebra::{Point3, Vector3};
use ply_rs::{parser, ply};
use std::io::BufRead;
use std::path::Path;
use std::str::FromStr;

pub type Real = f64;

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, PanOrbitCameraPlugin))
        .add_plugins(WireframePlugin)
        .add_systems(Startup, setup_camera_and_light)
        .add_systems(Startup, setup_scene)
        .run();
}

fn setup_scene(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let point_cloud = parse_file("./assets/xiaojiejie2_pcd.ply", true);

    let surface = poisson_reconstruction_impl::reconstruct_surface(&point_cloud);
    poisson_reconstruction_impl::spawn_mesh(&mut commands, &mut meshes, &mut materials, surface);

    println!();

    let surface = poissonrecon_impl::reconstruct_surface(&point_cloud);
    poissonrecon_impl::spawn_mesh(&mut commands, &mut meshes, &mut materials, surface);

    println!();
    println!("Done");
    println!();
    println!("mkazhdan's output is on the left, foresight mining's on the right.");
}

fn setup_camera_and_light(mut commands: Commands) {
    commands.spawn(PointLightBundle {
        point_light: PointLight {
            intensity: 100_000_000.0,
            range: 1000.,
            ..default()
        },
        transform: Transform::from_xyz(60.0, 50.0, 100.0),
        ..default()
    });
    commands.spawn((
        Camera3dBundle {
            transform: Transform::from_xyz(60.0, 25.0, 200.0)
                .looking_at(Vec3::new(60., 0., 0.), Vec3::Y),
            ..default()
        },
        PanOrbitCamera::default(),
    ));
}

#[derive(Default)]
struct VertexWithNormal {
    pos: Point3<Real>,
    normal: Vector3<Real>,
}

impl ply::PropertyAccess for VertexWithNormal {
    fn new() -> Self {
        Self::default()
    }

    fn set_property(&mut self, key: String, property: ply::Property) {
        match (key.as_ref(), property) {
            ("x", ply::Property::Float(v)) => self.pos.x = v as Real,
            ("y", ply::Property::Float(v)) => self.pos.y = v as Real,
            ("z", ply::Property::Float(v)) => self.pos.z = v as Real,
            ("nx", ply::Property::Float(v)) => self.normal.x = v as Real,
            ("ny", ply::Property::Float(v)) => self.normal.y = v as Real,
            ("nz", ply::Property::Float(v)) => self.normal.z = v as Real,
            _ => {}
        }
    }
}

fn parse_file(path: impl AsRef<Path>, ply: bool) -> Vec<VertexWithNormal> {
    let f = std::fs::File::open(path).unwrap();
    let mut f = std::io::BufReader::new(f);

    if ply {
        let vertex_parser = parser::Parser::<VertexWithNormal>::new();
        let header = vertex_parser.read_header(&mut f).unwrap();

        // Depending on the header, read the data into our structs..
        let mut vertex_list = Vec::new();
        for (_ignore_key, element) in &header.elements {
            // we could also just parse them in sequence, but the file format might change
            if element.name == "vertex" {
                vertex_list = vertex_parser
                    .read_payload_for_element(&mut f, element, &header)
                    .unwrap();
            }
        }
        vertex_list
    } else {
        let mut result = vec![];
        for line in f.lines().map_while(Result::ok) {
            let values: Vec<_> = line
                .split_whitespace()
                .map(|elt| f64::from_str(elt).unwrap())
                .collect();
            result.push(VertexWithNormal {
                pos: Point3::new(values[0], values[1], values[2]),
                normal: Vector3::new(values[3], values[4], values[5]),
            });
        }
        result
    }
}

mod poisson_reconstruction_impl {

    use super::VertexWithNormal;
    use bevy::pbr::wireframe::Wireframe;
    use bevy::prelude::*;
    use bevy::render::mesh::{Indices, PrimitiveTopology};
    use bevy::{asset::RenderAssetUsages, pbr::wireframe::WireframeColor};

    use poisson_reconstruction::PoissonReconstruction;
    use poisson_reconstruction::marching_cubes::MeshBuffers;

    pub fn reconstruct_surface(vertices: &[VertexWithNormal]) -> MeshBuffers {
        let points: Vec<_> = vertices.iter().map(|v| v.pos).collect();
        let normals: Vec<_> = vertices.iter().map(|v| v.normal).collect();

        println!("Running foresight mining's poisson implementation...");
        let start_time = std::time::Instant::now();
        let poisson =
            PoissonReconstruction::from_points_and_normals(&points, &normals, 0.0, 6, 6, 10);
        let buffers = poisson.reconstruct_mesh_buffers();
        let time_took = start_time.elapsed();
        println!(
            "foresight mining's poisson implementation took {:.3} seconds.",
            time_took.as_secs_f32()
        );
        buffers
    }

    pub fn spawn_mesh(
        commands: &mut Commands,
        meshes: &mut Assets<Mesh>,
        materials: &mut Assets<StandardMaterial>,
        points: MeshBuffers,
    ) {
        // Create the bevy mesh.
        let vertices: Vec<_> = points
            .vertices()
            .iter()
            .map(|pt| [pt.x as f32, pt.y as f32, pt.z as f32])
            .collect();
        let mut mesh = Mesh::new(
            PrimitiveTopology::TriangleList,
            RenderAssetUsages::default(),
        );
        mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, vertices);
        mesh.insert_indices(Indices::U32(points.indices().to_vec()));
        mesh.duplicate_vertices();
        mesh.compute_flat_normals();

        commands.spawn((
            Mesh3d(meshes.add(mesh)),
            MeshMaterial3d(materials.add(Color::srgb(1.0, 0.0, 1.0))),
            Transform::from_rotation(Quat::from_rotation_x(180.0f32.to_radians()))
                .with_translation([60.0, 0.0, 0.0].into()),
            Wireframe,
            WireframeColor {
                color: Color::srgb(0.3, 0.0, 0.3),
            },
        ));
    }
}

mod poissonrecon_impl {

    use super::VertexWithNormal;
    use bevy::asset::RenderAssetUsages;
    use bevy::pbr::wireframe::{Wireframe, WireframeColor};
    use bevy::prelude::*;
    use bevy::render::mesh::{Indices, PrimitiveTopology};
    use poissonrecon::{Mesh3D, PoissonParamsBuilder, SolutionParametersBuilder};

    pub fn reconstruct_surface(vertices: &[VertexWithNormal]) -> Mesh3D {
        let points: Vec<_> = vertices.iter().map(|v| v.pos).collect();
        let normals: Vec<_> = vertices.iter().map(|v| v.normal).collect();

        println!("Running mkazhdan's poisson implementation...");
        let start_time = std::time::Instant::now();
        let mesh = poissonrecon::reconstruct_surface(
            &points,
            &normals,
            &PoissonParamsBuilder::default()
                .solution_params(
                    SolutionParametersBuilder::default()
                        .point_weight(0.0)
                        .iters(10)
                        .depth(6)
                        .kernel_depth(6)
                        .build(),
                )
                .build(),
        )
        .unwrap();
        let time_took = start_time.elapsed();
        println!(
            "mkazhdan's poisson implementation took {:.3} seconds.",
            time_took.as_secs_f32()
        );
        mesh
    }

    pub fn spawn_mesh(
        commands: &mut Commands,
        meshes: &mut Assets<Mesh>,
        materials: &mut Assets<StandardMaterial>,
        surface: Mesh3D,
    ) {
        // Create the bevy mesh.
        let mut mesh = Mesh::new(
            PrimitiveTopology::TriangleList,
            RenderAssetUsages::default(),
        );
        mesh.insert_attribute(
            Mesh::ATTRIBUTE_POSITION,
            surface
                .vertices
                .iter()
                .map(|pt| [pt.x as f32, pt.y as f32, pt.z as f32])
                .collect::<Vec<_>>(),
        );
        mesh.insert_indices(Indices::U32(
            surface
                .triangles
                .iter()
                .flatten()
                .map(|&i| i as u32)
                .collect(),
        ));
        mesh.duplicate_vertices();
        mesh.compute_flat_normals();

        commands.spawn((
            Mesh3d(meshes.add(mesh)),
            MeshMaterial3d(materials.add(Color::srgb(1.0, 0.0, 1.0))),
            Transform::from_rotation(Quat::from_rotation_x(180.0f32.to_radians())),
            Wireframe,
            WireframeColor {
                color: Color::srgb(0.3, 0.0, 0.3),
            },
        ));
    }
}
