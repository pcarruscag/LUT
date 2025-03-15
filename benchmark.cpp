#include <fstream>
#include <iostream>
#include <omp.h>
#include <stdexcept>
#include <string>

#include "containers.hpp"
#include "query_map.hpp"

// Use the parallel (pstl + omp) build algorithm.
#define PARALLEL 1

#if PARALLEL
#include "build_map_par.hpp"
#include <tbb/global_control.h>
#else
#include "build_map_seq.hpp"
#endif

// Tolerance for tests.
constexpr RealT cTol = 1e-12;

int main() {
#if PARALLEL
  // Set TBB number of threads equal to OpenMP.
  tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism,
                                   omp_get_max_threads());
#else
  omp_set_num_threads(1);
#endif

  // Read the mesh (n tri, triangles, n pt, coords).
  std::ifstream f;
  f.open("mesh.su2");

  std::string line;
  getline(f, line); // NDIM
  getline(f, line);
  line = line.substr(6);

  const auto n_tri = static_cast<IntT>(std::stol(line));
  std::cout << n_tri << " triangles.\n";
  Matrix3i tris;
  tris.resize(n_tri, 3);

  IntT dummy;
  for (IntT i = 0; i < n_tri; ++i) {
    f >> dummy >> tris(i, 0) >> tris(i, 1) >> tris(i, 2) >> dummy;
  }

  getline(f, line);
  getline(f, line);
  line = line.substr(6);
  const auto n_pt = static_cast<IntT>(std::stol(line));
  std::cout << n_pt << " points.\n";
  VectorReal x, y;
  x.resize(n_pt);
  y.resize(n_pt);

  for (IntT i = 0; i < n_pt; ++i) {
    f >> x[i] >> y[i] >> dummy;
  }
  f.close();

  // Build the map.
  const auto t0 = omp_get_wtime();

  ReorderPoints<IntT>(tris, x, y);
  const auto t_order = omp_get_wtime();

  Matrix2i edge_pts, edge_faces;
  ExtractEdges<IntT>(tris, edge_pts, edge_faces);
  const auto t_edges = omp_get_wtime();

  TrapezoidalMap<VectorInt, VectorReal> map;
  BuildTrapezoidalMap<IntT>(edge_pts, x, y, map);
  const auto t_build = omp_get_wtime();

  std::cout << "Built map in " << (t_build - t0) << "s\n"
            << "  " << edge_pts.rows() << " edges and "
            << map.x_bands.size() - 1 << " bands\n"
            << "  Reorder in " << (t_order - t0) << "s\n"
            << "  Extract edges in " << (t_edges - t_order) << "s\n";

  // Manufacture a function.
  auto func = [](const auto x, const auto y) {
    return sin((x - y + 2) * 0.5 * M_PI);
  };

  VectorReal z;
  z.resize(n_pt);
#pragma omp parallel for
  for (IntT i = 0; i < n_pt; ++i) {
    z[i] = func(x[i], y[i]);
  }
  // Query the map at triangle centroids.
  VectorReal xi, yi, zi, zi_ref;
  xi.resize(n_tri);
  yi.resize(n_tri);
  zi.resize(n_tri);
  zi_ref.resize(n_tri);

#pragma omp parallel for
  for (IntT i = 0; i < n_tri; ++i) {
    const auto a = tris(i, 0);
    const auto b = tris(i, 1);
    const auto c = tris(i, 2);
    xi[i] = (x[a] + x[b] + x[c]) / 3;
    yi[i] = (y[a] + y[b] + y[c]) / 3;
    // Expected result of the interpolation.
    zi_ref[i] = (z[a] + z[b] + z[c]) / 3;

    // Test the parametric coordinates and containment test.
    const auto abc =
        TriangleCoords(x[a], y[a], x[b], y[b], x[c], y[c], xi[i], yi[i]);
    if (!InTriangle(abc)) {
      throw std::runtime_error("Point was expected to be inside the triangle.");
    }
    constexpr auto third = RealT(1) / 3;
    if (fabs(abc[0] - third) > cTol || fabs(abc[1] - third) > cTol) {
      throw std::runtime_error("Wrong parametric coordinates.");
    }
  }

  auto t_query = -omp_get_wtime();
#pragma omp parallel for
  for (IntT i = 0; i < n_tri; ++i) {
    zi[i] = Interpolate<IntT>(map, edge_faces, tris, xi[i], yi[i], x, y, z)[0];
  }
  t_query += omp_get_wtime();
  std::cout << "Queried map in " << t_query << "s\n";

  RealT max_err = 0, avg_err = 0;
  for (IntT i = 0; i < n_tri; ++i) {
    const auto val = func(xi[i], yi[i]);
    const auto err = fabs(val - zi[i]);
    max_err = fmax(max_err, err);
    avg_err += err;
    if (fabs(zi[i] - zi_ref[i]) > cTol) {
      throw std::runtime_error(
          "Wrong interpolated value (incorrect triangle).");
    }
  }
  avg_err /= static_cast<RealT>(n_tri);
  std::cout << "  Max/Avg error: " << max_err << "/" << avg_err << "\n";

  return 0;
}