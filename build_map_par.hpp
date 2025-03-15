#include <algorithm>
#include <array>
#include <execution>
#include <iostream>
#include <iterator>
#include <numeric>
#include <omp.h>
#include <tuple>
#include <utility>
#include <vector>

template <typename Vector> struct ScalarType;

template <typename T> struct ScalarType<std::vector<T>> {
  using type = T;
};

template <typename Vector>
using ScalarTypeT = typename ScalarType<Vector>::type;

/// The map is defined by the limits of the bands in the x direction
/// and a CSR of the edge IDs in each band, sorted by the midpoint of
/// the edge segment in that band.
template <typename VectorInt, typename VectorReal> struct TrapezoidalMap {
  VectorInt offsets, edge_id;
  VectorReal x_bands, edge_y;
};

/// Orders points in order of ascending x coordinates and updates the definition
/// of triangles.
template <typename Int, typename Matrix3, typename VectorReal>
void ReorderPoints(Matrix3 &triangles, VectorReal &x, VectorReal &y) {
  const auto n_pts = static_cast<Int>(x.size());

  // Sort order of x values.
  std::vector<Int> perm(n_pts);
#pragma omp parallel for
  for (Int i = 0; i < n_pts; ++i) {
    perm[i] = i;
  }
  std::sort(std::execution::par, perm.begin(), perm.end(),
            [&x, &y](const auto i, const auto j) {
              return x[i] != x[j] ? x[i] < x[j] : y[i] < y[j];
            });

  auto reorder = [n_pts, &perm](const auto &v) {
    VectorReal tmp(n_pts);
#pragma omp parallel for
    for (Int i = 0; i < n_pts; ++i) {
      tmp[i] = v[perm[i]];
    }
    return tmp;
  };
  x = reorder(x);
  y = reorder(y);

  std::vector<Int> inv_perm(n_pts);
#pragma omp parallel
  {
#pragma omp for
    for (Int i = 0; i < n_pts; ++i) {
      inv_perm[perm[i]] = i;
    }
#pragma omp for
    for (Int i = 0; i < static_cast<Int>(triangles.rows()); ++i) {
      for (Int j = 0; j < 3; ++j) {
        triangles(i, j) = inv_perm[triangles(i, j)];
      }
    }
  }
}

/// Extracts unique edges from a set of triangles, edges are defined by two
/// point IDs and can be shared by up to two faces (triangles), boundary edges
/// only have one adjacent face (the second face ID will be < 0).
template <typename Int, typename Matrix3, typename Matrix2>
void ExtractEdges(const Matrix3 &triangles, Matrix2 &edge_pts,
                  Matrix2 &edge_faces) {
  // Extract edges from the triangles.
  std::vector<std::array<Int, 3>> edges;
  edges.resize(3 * triangles.rows());

#pragma omp parallel for
  for (Int i_tri = 0; i_tri < static_cast<Int>(triangles.rows()); ++i_tri) {
    for (Int i = 0; i < 3; ++i) {
      const auto j = (i + 1) % 3;
      const auto i_pt = std::min(triangles(i_tri, i), triangles(i_tri, j));
      const auto j_pt = std::max(triangles(i_tri, i), triangles(i_tri, j));
      edges[3 * i_tri + i] = {i_pt, j_pt, i_tri};
    }
  }
  // Sort to identify duplicates.
  std::sort(std::execution::par, edges.begin(), edges.end(),
            [](const auto &a, const auto &b) {
              return a[0] != b[0] ? (a[0] < b[0]) : (a[1] < b[1]);
            });

  // Mark unique edges and determine where the unique edges should be stored.
  std::vector<Int> new_edge_offset(edges.size());
  new_edge_offset[0] = 0;

#pragma omp parallel for
  for (Int i = 1; i < static_cast<Int>(edges.size()); ++i) {
    auto is_new = [](const auto &a, const auto &b) {
      return static_cast<Int>(a[0] != b[0] || a[1] != b[1]);
    };
    new_edge_offset[i] = is_new(edges[i], edges[i - 1]);
  }
  std::inclusive_scan(std::execution::par, new_edge_offset.begin(),
                      new_edge_offset.end(), new_edge_offset.begin());
  // Add 1 because we set new_edge_offset[0] to 0.
  const auto n_edges = new_edge_offset.back() + 1;

  // Map edge points and edge faces.
  edge_pts.resize(n_edges, 2);
  edge_faces.resize(n_edges, 2);

  auto new_edge = [&](const Int i, const auto &edge) {
    edge_pts(i, 0) = edge[0];
    edge_pts(i, 1) = edge[1];
    edge_faces(i, 0) = edge[2];
  };
  new_edge(0, edges[0]);

#pragma omp parallel
  {
    // Initialize the second face of the edge to avoid race conditions.
    // One thread may call new_edge while other tries to set the second face.
#pragma omp for
    for (Int i = 0; i < n_edges; ++i) {
      edge_faces(i, 1) = -1;
    }
#pragma omp for
    for (Int i = 1; i < static_cast<Int>(edges.size()); ++i) {
      const auto pos = new_edge_offset[i];
      if (pos != new_edge_offset[i - 1]) {
        new_edge(pos, edges[i]);
      } else {
        // Map the second face.
        edge_faces(pos, 1) = edges[i][2];
      }
    }
  }
}

/// Given the x coordinates of all points, detects the number of bands
/// (unique x values - 1). Returns the number of bands, their limits,
/// and a map from point ID to band ID.
template <typename Int, typename Vector> auto DetectBands(const Vector &x) {
  // Since there could be duplicate x values, the number of bands is not n-1.
  std::vector<Int> band(x.size());
  band[0] = 0;
#pragma omp parallel for
  for (Int i = 1; i < static_cast<Int>(x.size()); ++i) {
    band[i] = static_cast<Int>(x[i] != x[i - 1]);
  }
  // Determine where unique x values should be stored to define the bands.
  std::inclusive_scan(std::execution::par, band.begin(), band.end(),
                      band.begin());
  // Do not add 1 despite band[0] == 0 because the number of bands is the
  // number of unique x values minus 1.
  const auto n_bands = band.back();

  // The inverse permutation maps point indices (including duplicates) to the
  // band index.
  std::vector<Int> pt_to_band(x.size());
  Vector x_bands(n_bands + 1);

#pragma omp parallel for
  for (Int i = 0; i < static_cast<Int>(x.size()); ++i) {
    const auto pos = band[i];
    if (i == 0 || pos != band[i - 1]) {
      x_bands[pos] = x[i];
    }
    pt_to_band[i] = pos;
  }
  return std::make_tuple(n_bands, std::move(pt_to_band), std::move(x_bands));
}

/// Builds a trapezoidal map for a set of edges.
template <typename Int, typename Matrix2, typename VectorReal,
          typename VectorInt>
void BuildTrapezoidalMap(const Matrix2 &edge_pts, const VectorReal &x,
                         const VectorReal &y,
                         TrapezoidalMap<VectorInt, VectorReal> &map) {
  auto &x_bands = map.x_bands;
  auto &offsets = map.offsets;
  auto &edge_id = map.edge_id;
  auto &edge_y = map.edge_y;

  const auto [n_bands, pt_to_band, bands] = DetectBands<Int>(x);
  x_bands = std::move(bands);

  // Count the number of edges per band.
  auto &counts = offsets;
  counts.clear();
  counts.resize(n_bands + 1, 0);
#pragma omp parallel for
  for (Int i = 0; i < static_cast<Int>(edge_pts.rows()); ++i) {
    const auto band_0 = pt_to_band[edge_pts(i, 0)];
    const auto band_1 = pt_to_band[edge_pts(i, 1)];
    for (auto j = std::min(band_0, band_1); j < std::max(band_0, band_1); ++j) {
#pragma omp atomic
      ++counts[j + 1];
    }
  }
  // Convert the counts to offsets.
  std::inclusive_scan(std::execution::par, offsets.begin(), offsets.end(),
                      offsets.begin());

  auto t2 = -omp_get_wtime();

  // For each band store the edge ids that cross it, and the mid-point y
  // coordinate of the segment.
  edge_id.resize(offsets.back());
  edge_y.resize(offsets.back());
  auto pos = offsets;
#pragma omp parallel
  {
    // const auto thread = omp_get_thread_num();
    // const auto bands_per_thread =
    //     (n_bands + omp_get_num_threads() - 1) / omp_get_num_threads();
    // const auto band_begin = thread * bands_per_thread;
    // const auto band_end = std::min(band_begin + bands_per_thread, n_bands);

    // auto process_band = [&](const auto band) {
    //   return band >= band_begin && band < band_end;
    // };
#pragma omp for
    for (Int i_edge = 0; i_edge < static_cast<Int>(edge_pts.rows()); ++i_edge) {
      const auto pt_0 = edge_pts(i_edge, 0);
      const auto pt_1 = edge_pts(i_edge, 1);
      const auto band_0 = pt_to_band[pt_0];
      const auto band_1 = pt_to_band[pt_1];
      // Vertical edges are not mapped.
      if (band_0 == band_1)
        continue;
      const auto x_0 = x[pt_0];
      const auto y_0 = y[pt_0];
      const auto dy_dx = (y[pt_1] - y_0) / (x[pt_1] - x_0);
      for (auto j = std::min(band_0, band_1); j < std::max(band_0, band_1);
           ++j) {
        // if (!process_band(j))
        //   continue;
        Int pj{};
#pragma omp atomic capture
        pj = pos[j]++;
        edge_id[pj] = i_edge;
        const auto x_mid = (x_bands[j] + x_bands[j + 1]) / 2;
        edge_y[pj] = y_0 + dy_dx * (x_mid - x_0);
      }
    }
  }
  decltype(pos)().swap(pos);

  t2 += omp_get_wtime();
  std::cout << t2 << std::endl;

  auto t3 = -omp_get_wtime();

  // Sort the edges in each band by y coordinate.
#pragma omp parallel
  {
    std::vector<std::pair<Int, ScalarTypeT<VectorReal>>> tmp;
#pragma omp for
    for (Int i = 0; i < n_bands; ++i) {
      const auto begin = offsets[i];
      const auto end = offsets[i + 1];
      tmp.resize(end - begin);
      for (auto k = begin; k < end; ++k) {
        tmp[k - begin] = {edge_id[k], edge_y[k]};
      }
      std::sort(tmp.begin(), tmp.end(),
                [](const auto a, const auto b) { return a.second < b.second; });
      for (auto k = begin; k < end; ++k) {
        edge_id[k] = tmp[k - begin].first;
        edge_y[k] = tmp[k - begin].second;
      }
    }
  }

  t3 += omp_get_wtime();
  std::cout << t3 << std::endl;
}
