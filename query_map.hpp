#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <utility>

/// Returns the edge IDs below and above a query point. If the point is
/// below/above the lowest/highest edge the first/second ID will be < 0.
template <typename Int, typename Map, typename Real>
auto QueryTrapezoidalMap(const Map &map, const Real &x, const Real &y) {
  // Find the band.
  const auto i = [&]() {
    const auto &x_bands = map.x_bands;
    const auto it = std::lower_bound(x_bands.begin(), x_bands.end(), x);
    const auto d = static_cast<Int>(std::distance(x_bands.begin(), it));
    const auto n = static_cast<Int>(x_bands.size() - 2);
    return std::min(std::max<Int>(0, d - 1), n);
  }();

  // Find the edge above.
  const auto begin = map.offsets[i];
  const auto end = map.offsets[i + 1];
  const auto it =
      std::lower_bound(map.edge_y.begin() + begin, map.edge_y.begin() + end, y);
  const auto j = std::distance(map.edge_y.begin(), it);
  if (j == begin) {
    return std::make_pair(Int{-1}, map.edge_id[j]);
  }
  if (j == end) {
    return std::make_pair(map.edge_id[j - 1], Int{-1});
  }
  return std::make_pair(map.edge_id[j - 1], map.edge_id[j]);
}

/// Returns the ids of the triangles adjacent to 2 query edges.
/// Up to 3 trianges ids are returned.
template <typename Int, typename Matrix2>
auto AdjacentTriangles(const Int edge_0, const Int edge_1,
                       const Matrix2 &edge_faces) {
  std::array<Int, 3> tris = {-1, -1, -1};
  Int pos = 0;
  auto insert = [&tris, &pos](const Int t) {
    if (t < 0)
      return;
    for (Int i = 0; i < pos; ++i) {
      if (t == tris[i])
        return;
    }
    tris[pos++] = t;
  };
  auto get_tris = [&edge_faces](const Int e) {
    if (e < 0)
      return std::array{Int{-1}, Int{-1}};
    return std::array{edge_faces(e, 0), edge_faces(e, 1)};
  };
  for (const auto e : {edge_0, edge_1}) {
    for (const auto t : get_tris(e))
      insert(t);
  }
  return tris;
}

/// Computes a,b,c such that x = a x0 + b x1 + c x2.
template <typename Real>
std::array<Real, 3>
TriangleCoords(const Real &x0, const Real &y0, const Real &x1, const Real &y1,
               const Real &x2, const Real &y2, const Real &x, const Real &y) {
  const auto dx1 = x1 - x0;
  const auto dy1 = y1 - y0;
  const auto dx2 = x2 - x0;
  const auto dy2 = y2 - y0;

  auto cross = [](const Real &ux, const Real &uy, const Real &vx,
                  const Real &vy) { return ux * vy - uy * vx; };
  const auto inv_det = 1 / cross(dx1, dy1, dx2, dy2);
  const auto a = (cross(x, y, dx2, dy2) - cross(x0, y0, dx2, dy2)) * inv_det;
  const auto b = (cross(x0, y0, dx1, dy1) - cross(x, y, dx1, dy1)) * inv_det;
  return {1 - a - b, a, b};
}

/// Returns true if a point is inside a tringle based on the result from
/// TriangleCoords.
template <typename Real>
bool InTriangle(const std::array<Real, 3> &tri_coords) {
  return fmin(fmin(tri_coords[0], tri_coords[1]), tri_coords[2]) >= 0;
}

/// Interpolate the value(s) of zs at point (x, y).
template <typename Int, typename Map, typename Matrix2, typename Matrix3,
          typename Real, typename VectorReal, typename... VectorReals>
auto Interpolate(const Map &map, const Matrix2 &edge_faces,
                 const Matrix3 &triangles, const Real &x, const Real &y,
                 const VectorReal &xs, const VectorReal &ys,
                 const VectorReals &...zs) {
  const auto [e_0, e_1] = QueryTrapezoidalMap<Int>(map, x, y);
  const auto q_tris = AdjacentTriangles<Int>(e_0, e_1, edge_faces);
  std::array<Real, sizeof...(zs)> z{};
  for (const auto t : q_tris) {
    if (t < 0)
      continue;
    auto coords = [&](const Int p) { return std::make_pair(xs[p], ys[p]); };
    const auto p0 = triangles(t, 0);
    const auto p1 = triangles(t, 1);
    const auto p2 = triangles(t, 2);
    const auto [x0, y0] = coords(p0);
    const auto [x1, y1] = coords(p1);
    const auto [x2, y2] = coords(p2);
    const auto abc = TriangleCoords(x0, y0, x1, y1, x2, y2, x, y);
    if (InTriangle(abc)) {
      z = {abc[0] * zs[p0] + abc[1] * zs[p1] + abc[2] * zs[p2]...};
      break;
    }
  }
  return z;
}