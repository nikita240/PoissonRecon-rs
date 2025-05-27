#include "poisson_wrapper.h"
#include <cstdlib>
#include <vector>

using namespace PoissonRecon;

// Custom stream that reads from memory instead of from disk
class MemoryOrientedSampleStream
    : public Reconstructor::InputOrientedSampleStream<Real, 3> {
public:
  MemoryOrientedSampleStream(const PointNormalData *data) : _data(data) {}

  void reset() override { _current = 0; }

  bool read(Point<Real, 3> &p, Point<Real, 3> &n) override {
    if (_current >= _data->count)
      return false;

    p = _data->points[_current];
    n = _data->normals[_current];

    _current++;
    return true;
  }

private:
  const PointNormalData *_data;
  size_t _current;
};

// Output stream for vertices
class VectorVertexStream
    : public Reconstructor::OutputLevelSetVertexStream<Real, 3> {
public:
  VectorVertexStream(std::vector<Real> *vertices) : _vertices(vertices) {}

  size_t write(const Reconstructor::Position<Real, 3> &p,
               const Reconstructor::Gradient<Real, 3> &,
               const Reconstructor::Weight<Real> &) override {
    size_t idx = size();
    // NOTE: We do not currently support vertex attributes, so we ignore the
    // gradient and weight.
    for (int i = 0; i < 3; i++) {
      _vertices->push_back(p[i]);
    }
    return idx;
  }

  size_t size() const override { return _vertices->size() / 3; }

private:
  std::vector<Real> *_vertices;
};

// Output stream for faces
class VectorFaceStream : public Reconstructor::OutputFaceStream<2> {
public:
  VectorFaceStream(std::vector<int> *indices) : _indices(indices) {}

  size_t write(const std::vector<node_index_type> &polygon) override {
    size_t idx = size();
    // NOTE: This assumes that the polygon is a triangle (3 vertices).
    // The Poisson bindings currently do not support quads and other polygon
    // types.
    for (int i = 0; i < 3; i++) {
      _indices->push_back(polygon[i]);
    }
    return idx;
  }

  size_t size() const override { return _indices->size() / 3; }

private:
  std::vector<int> *_indices;
};

// The main reconstruction function
template <unsigned int Degree, BoundaryType BType>
MeshData *poisson_reconstruct(const PointNormalData *data,
                              CPoissonParams params) {
  // Create sample stream from the provided points and normals
  MemoryOrientedSampleStream sampleStream(data);

  // Configure parameters
  Reconstructor::Poisson::SolutionParameters<Real> solverParams =
      params.solution_params;
  Reconstructor::LevelSetExtractionParameters extractionParams =
      params.level_set_extraction_params;

  // The tensor-product finite-elements signatures
  static const unsigned int FEMSig =
      FEMDegreeAndBType<Degree, BType>::Signature;
  using FEMSigs = IsotropicUIntPack<3, FEMSig>;
  using Implicit = Reconstructor::Implicit<Real, 3, FEMSigs>;
  using Solver = Reconstructor::Poisson::Solver<Real, 3, FEMSigs>;

  // Run the solver
  Implicit *implicit;
  try {
    implicit = Solver::Solve(sampleStream, solverParams);
  } catch (const std::exception &e) {
    std::cerr << "Error during reconstruction: " << e.what() << std::endl;
    return nullptr;
  }
  if (!implicit)
    return nullptr;

  // Set up output vectors
  std::vector<Real> *vertices = new std::vector<Real>;
  std::vector<int> *indices = new std::vector<int>;

  // Create streams to capture the output
  VectorVertexStream vertexStream(vertices);
  VectorFaceStream faceStream(indices);

  // Extract level set (the actual surface)
  try {
    implicit->extractLevelSet(vertexStream, faceStream, extractionParams);
  } catch (const std::exception &e) {
    std::cerr << "Error during level set extraction: " << e.what() << std::endl;
    delete implicit;
    return nullptr;
  }

  // Clean up
  delete implicit;

  // We are intentionally leaking the data.
  MeshData *mesh = new MeshData();
  mesh->vertices = vertices->data();
  mesh->indices = indices->data();
  mesh->vertex_count = vertices->size() / 3;
  mesh->index_count = indices->size() / 3;

  // Delete the vectors without deleting the data.
  // It's up to the calling Rust code to delete the data after it has copied it
  // out.
  free(vertices);
  free(indices);

  return mesh;
}

template <unsigned int Degree>
MeshData *poisson_reconstruct(const PointNormalData *data,
                              CPoissonParams params,
                              CBoundaryType boundary_type) {
#ifdef EN_BOUNDARY_FREE
  if (boundary_type == CBoundaryType::Free) {
    return poisson_reconstruct<Degree, BOUNDARY_FREE>(data, params);
  }
#endif
#ifdef EN_BOUNDARY_DIRICHLET
  if (boundary_type == CBoundaryType::Dirichlet) {
    return poisson_reconstruct<Degree, BOUNDARY_DIRICHLET>(data, params);
  }
#endif
#ifdef EN_BOUNDARY_NEUMANN
  if (boundary_type == CBoundaryType::Neumann) {
    return poisson_reconstruct<Degree, BOUNDARY_NEUMANN>(data, params);
  }
#endif

  return nullptr;
}

// Our "unrolling" entrypoint
extern "C" MeshData *poisson_reconstruct(const PointNormalData *data,
                                         CPoissonParams params, CDegree degree,
                                         CBoundaryType boundary_type) {
#ifdef EN_DEGREE_1
  if ((uint8_t)degree == 1) {
    return poisson_reconstruct<1>(data, params, boundary_type);
  }
#endif
#ifdef EN_DEGREE_2
  if ((uint8_t)degree == 2) {
    return poisson_reconstruct<2>(data, params, boundary_type);
  }
#endif
#ifdef EN_DEGREE_3
  if ((uint8_t)degree == 3) {
    return poisson_reconstruct<3>(data, params, boundary_type);
  }
#endif

  return nullptr;
}

extern "C" void free_mesh_data(MeshData *mesh) {
  if (mesh) {
    free(mesh->vertices);
    free(mesh->indices);
    free(mesh);
  }
}

extern "C" CPoissonParams get_default_poisson_params() {
  CPoissonParams params;
  return params;
}
