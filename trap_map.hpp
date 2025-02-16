#include <iterator>
#include <vector>
#include <utility>
#include <array>
#include <algorithm>
#include <numeric>
#include <tuple>

template <typename Vector>
struct ScalarType;

template <typename T>
struct ScalarType<std::vector<T>> {
    using type = T;
};

template <typename Vector>
using ScalarTypeT = typename ScalarType<Vector>::type;

/// The map is defined by the limits of the bands in the x direction
/// and a CSR of the edge IDs in each band, sorted by the midpoint of
/// the edge segment in that band.
template <typename VectorInt, typename VectorReal>
struct TrapezoidalMap {
    VectorInt offsets, edge_id;
    VectorReal x_bands, edge_y;
};

/// Extracts unique edges from a set of triangles, edges are defined by two
/// point IDs and can be shared by up to two faces (triangles), boundary edges
/// only have one adjacent face (the second face ID will be < 0).
template <typename Int, typename Matrix3, typename Matrix2>
void ExtractEdges(const Matrix3& triangles, Matrix2& edge_pts, Matrix2& edge_faces) {
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
    std::sort(edges.begin(), edges.end(), [](const auto& a, const auto& b) {
        return a[0] != b[0] ? (a[0] < b[0]) : (a[1] < b[1]);
    });
    // Count unique edges.
    Int n_edges = 1;
    auto is_equal = [](const auto& a, const auto& b) {
        return a[0] == b[0] && a[1] == b[1];
    };
    for (Int i = 1; i < static_cast<Int>(edges.size()); ++i) {
        n_edges += static_cast<Int>(!is_equal(edges[i], edges[i - 1]));
    }
    // Map edge points and edge faces.
    edge_pts.resize(n_edges, 2);
    edge_faces.resize(n_edges, 2);
    Int pos = 0;
    auto new_edge = [&](const auto& edge) {
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
template <typename Int, typename Vector>
auto DetectBands(const Vector& x) {
    // Sort order of x values.
    std::vector<Int> perm(x.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&x](const auto i, const auto j) {
        return x[i] < x[j];
    });
    // Since there could be duplicate x values, the number of bands is not n-1.
    std::vector<Int> band(perm.size());
    band[0] = 0;
    Int n_bands = 0;
    for (Int i = 1; i < static_cast<Int>(perm.size()); ++i) {
        n_bands += static_cast<Int>(x[perm[i]] != x[perm[i - 1]]);
        band[i] = n_bands;
    }
    Vector x_bands(n_bands + 1);
    Int pos = 0;
    x_bands[pos] = x[perm[0]];
    for (Int i = 1; i < static_cast<Int>(perm.size()); ++i) {
        if (x[perm[i]] != x_bands[pos]) {
            x_bands[++pos] = x[perm[i]];
        }
    }
    // The inverse permutation maps point indices to the band index.
    std::vector<Int> pt_to_band(perm.size());
    for (Int i = 0; i < static_cast<Int>(perm.size()); ++i) {
        pt_to_band[perm[i]] = band[i];
    }
    return std::make_tuple(n_bands, std::move(pt_to_band), std::move(x_bands));
}

/// Builds a trapezoidal map for a set of edges.
template <typename Int, typename Matrix2, typename VectorReal, typename VectorInt>
void BuildTrapezoidalMap(const Matrix2& edge_pts, const VectorReal& x, const VectorReal& y,
                         TrapezoidalMap<VectorInt, VectorReal>& map) {
    auto& x_bands = map.x_bands;
    auto& offsets = map.offsets;
    auto& edge_id = map.edge_id;
    auto& edge_y = map.edge_y;

    const auto [n_bands, pt_to_band, bands] = DetectBands<Int>(x);
    x_bands = std::move(bands);

    // Count the number of edges per band.
    auto& counts = offsets;
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
        if (band_0 == band_1) continue;
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
        std::sort(tmp.begin(), tmp.end(), [](const auto a, const auto b) {
            return a.second < b.second;
        });
        for (auto k = begin; k < end; ++k) {
            edge_id[k] = tmp[k - begin].first;
            edge_y[k] = tmp[k - begin].second;
        }
    }
}

/// Returns the edge IDs below and above a query point. If the point is
/// below/above the lowest/highest edge the first/second ID will be < 0.
template <typename Int, typename Map, typename Real>
auto QueryTrapezoidalMap(const Map& map, const Real& x, const Real& y) {
    // Find the band.
    const auto i = [&]() {
        const auto& x_bands = map.x_bands;
        const auto it = std::lower_bound(x_bands.begin(), x_bands.end(), x);
        const auto d = static_cast<Int>(std::distance(x_bands.begin(), it));
        const auto n = static_cast<Int>(x_bands.size() - 2);
        return std::min(std::max<Int>(0, d - 1), n);
    }();

    // Find the edge above.
    const auto begin = map.offsets[i];
    const auto end = map.offsets[i + 1];
    const auto it = std::lower_bound(map.edge_y.begin() + begin,
                                     map.edge_y.begin() + end, y);
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
                       const Matrix2& edge_faces) {
    std::array<Int, 3> tris = {-1, -1, -1};
    Int pos = 0;
    auto insert = [&tris, &pos](const Int t) {
        if (t < 0) return;
        for (Int i = 0; i < pos; ++i) {
            if (t == tris[i]) return;
        }
        tris[pos++] = t;
    };
    auto get_tris = [&edge_faces](const Int e) {
        if (e < 0) return std::array{Int{-1}, Int{-1}};
        return std::array{edge_faces(e, 0), edge_faces(e, 1)};
    };
    for (const auto e : {edge_0, edge_1}) {
        for (const auto t : get_tris(e)) insert(t);
    }
    return tris;
}

/// Computes a,b,c such that x = a x0 + b x1 + c x2.
template <typename Real>
std::array<Real, 3> TriangleCoords(
        const Real& x0, const Real& y0, const Real& x1, const Real& y1,
        const Real& x2, const Real& y2, const Real& x, const Real& y) {
    const auto dx1 = x1 - x0;
    const auto dy1 = y1 - y0;
    const auto dx2 = x2 - x0;
    const auto dy2 = y2 - y0;

    auto cross = [](const Real& ux, const Real& uy,
                    const Real& vx, const Real& vy) {
        return ux * vy - uy * vx;
    };
    const auto inv_det = 1 / cross(dx1, dy1, dx2, dy2);
    const auto a = (cross(x, y, dx2, dy2) - cross(x0, y0, dx2, dy2)) * inv_det;
    const auto b = (cross(x0, y0, dx1, dy1) - cross(x, y, dx1, dy1)) * inv_det;
    return {1 - a - b, a, b};
}

/// Returns true if a point is inside a tringle based on the result from
/// TriangleCoords.
template <typename Real>
bool InTriangle(const std::array<Real, 3>& tri_coords) {
    return tri_coords[0] <= 1 && tri_coords[1] >= 0 && tri_coords[2] >= 0;
}

/// Interpolate the value(s) of zs at point (x, y).
template <typename Int, typename Map, typename Matrix2, typename Matrix3,
          typename Real, typename VectorReal,  typename... VectorReals>
auto Interpolate(const Map& map, const Matrix2& edge_faces,
                 const Matrix3& triangles, const Real& x, const Real& y,
                 const VectorReal& xs, const VectorReal& ys,
                 const VectorReals&... zs) {
    const auto [e_0, e_1] = QueryTrapezoidalMap<Int>(map, x, y);
    const auto q_tris = AdjacentTriangles<Int>(e_0, e_1, edge_faces);
    std::array<Real, sizeof...(zs)> z{};
    for (const auto t : q_tris) {
        if (t < 0) continue;
        auto coords = [&](const Int p) {
            return std::make_pair(xs[p], ys[p]);
        };
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