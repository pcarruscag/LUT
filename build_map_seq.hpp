#include <algorithm>
#include <array>
#include <iterator>
#include <numeric>
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
  std::iota(perm.begin(), perm.end(), 0);
  std::sort(perm.begin(), perm.end(), [&x, &y](const auto i, const auto j) {
    return x[i] != x[j] ? x[i] < x[j] : y[i] < y[j];
  });

  auto reorder = [n_pts, &perm](const auto &v) {
    VectorReal tmp(n_pts);
    for (Int i = 0; i < n_pts; ++i) {
      tmp[i] = v[perm[i]];
    }
    return tmp;
  };
  x = reorder(x);
  y = reorder(y);

  std::vector<Int> inv_perm(n_pts);
  for (Int i = 0; i < n_pts; ++i) {
    inv_perm[perm[i]] = i;
  }
  for (Int i = 0; i < static_cast<Int>(triangles.rows()); ++i) {
    for (Int j = 0; j < 3; ++j) {
      triangles(i, j) = inv_perm[triangles(i, j)];
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
  for (Int i_tri = 0; i_tri < static_cast<Int>(triangles.rows()); ++i_tri) {
    for (Int i = 0; i < 3; ++i) {
      const auto j = (i + 1) % 3;
      const auto i_pt = std::min(triangles(i_tri, i), triangles(i_tri, j));
      const auto j_pt = std::max(triangles(i_tri, i), triangles(i_tri, j));
      edges[3 * i_tri + i] = {i_pt, j_pt, i_tri};
    }
  }
  // Sort to identify duplicates.
  std::sort(edges.begin(), edges.end(), [](const auto &a, const auto &b) {
    return a[0] != b[0] ? (a[0] < b[0]) : (a[1] < b[1]);
  });
  // Count unique edges.
  Int n_edges = 1;
  auto is_equal = [](const auto &a, const auto &b) {
    return a[0] == b[0] && a[1] == b[1];
  };
  for (Int i = 1; i < static_cast<Int>(edges.size()); ++i) {
    n_edges += static_cast<Int>(!is_equal(edges[i], edges[i - 1]));
  }
  // Map edge points and edge faces.
  edge_pts.resize(n_edges, 2);
  edge_faces.resize(n_edges, 2);
  Int pos = 0;
  auto new_edge = [&](const auto &edge) {
    edge_pts(pos, 0) = edge[0];
    edge_pts(pos, 1) = edge[1];
    edge_faces(pos, 0) = edge[2];
    edge_faces(pos, 1) = -1;
    ++pos;
  };
  new_edge(edges[0]);
  for (Int i = 1; i < static_cast<Int>(edges.size()); ++i) {
    if (is_equal(edges[i], edges[i - 1])) {
      edge_faces(pos - 1, 1) = edges[i][2];
    } else {
      new_edge(edges[i]);
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
  Int n_bands = 0;
  for (Int i = 1; i < static_cast<Int>(x.size()); ++i) {
    n_bands += static_cast<Int>(x[i] != x[i - 1]);
    band[i] = n_bands;
  }
  Vector x_bands(n_bands + 1);
  Int pos = 0;
  x_bands[pos] = x[0];
  for (Int i = 1; i < static_cast<Int>(x.size()); ++i) {
    if (x[i] != x_bands[pos]) {
      x_bands[++pos] = x[i];
    }
  }
  return std::make_tuple(n_bands, std::move(band), std::move(x_bands));
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
  for (Int i = 0; i < static_cast<Int>(edge_pts.rows()); ++i) {
    const auto band_0 = pt_to_band[edge_pts(i, 0)];
    const auto band_1 = pt_to_band[edge_pts(i, 1)];
    for (auto j = std::min(band_0, band_1); j < std::max(band_0, band_1); ++j) {
      ++counts[j + 1];
    }
  }
  // Convert the counts to offsets.
  for (Int i = 2; i < static_cast<Int>(offsets.size()); ++i) {
    offsets[i] += offsets[i - 1];
  }
  // For each band store the edge ids that cross it,
  // and the mid-point y coordinate of the segment.
  edge_id.resize(offsets.back());
  edge_y.resize(offsets.back());
  auto pos = offsets;
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
    for (auto j = std::min(band_0, band_1); j < std::max(band_0, band_1); ++j) {
      edge_id[pos[j]] = i_edge;
      const auto x_mid = (x_bands[j] + x_bands[j + 1]) / 2;
      edge_y[pos[j]] = y_0 + dy_dx * (x_mid - x_0);
      ++pos[j];
    }
  }
  decltype(pos)().swap(pos);
  // Sort the edges in each band by y coordinate.
  std::vector<std::pair<Int, ScalarTypeT<VectorReal>>> tmp;
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
