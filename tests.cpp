#include <iostream>

#include "containers.hpp"
#include "trap_map.hpp"


int main() {
    /*
     * 0 -- 1
     * | \  |
     * |  \ |
     * 2 -- 3
     */
    Matrix3i tris;
    tris.resize(2, 3);
    tris(0, 0) = 0;
    tris(0, 1) = 3;
    tris(0, 2) = 1;
    tris(1, 0) = 0;
    tris(1, 1) = 2;
    tris(1, 2) = 3;

    const VectorReal x = {0, 1, -0.5, 1};
    const VectorReal y = {0.5, 0.5, -0.5, -0.5};

    Matrix2i edge_pts, edge_faces;
    ExtractEdges<IntT>(tris, edge_pts, edge_faces);

    std::cout << "\nPoints\n";
    for (size_t i = 0; i < edge_pts.rows(); ++i) {
        std::cout << edge_pts(i, 0) << "  " << edge_pts(i, 1) << '\n';
    }
    std::cout << "\nFaces\n";
    for (size_t i = 0; i < edge_faces.rows(); ++i) {
        std::cout << edge_faces(i, 0) << "  " << edge_faces(i, 1) << '\n';
    }

    TrapezoidalMap<VectorInt, VectorReal> map;
    BuildTrapezoidalMap<IntT>(edge_pts, x, y, map);

    std::cout << "\nBand Limits\n";
    for (const auto& x : map.x_bands) {
        std::cout << x << '\n';
    }
    std::cout << "\nBand Edges\n";
    for (size_t i = 0; i < map.offsets.size() - 1; ++i) {
        std::cout << i << '\n';
        for (auto j = map.offsets[i]; j < map.offsets[i + 1]; ++j) {
            std::cout << map.edge_id[j] << "  " << map.edge_y[j] << '\n';
        }
    }

    const auto px = 0.25;
    const auto py = 0.249;

    std::cout << "\nQuery\n";
    const auto [e_0, e_1] = QueryTrapezoidalMap<IntT>(map, px, py);
    std::cout << e_0 << "  " << e_1 << '\n';

    std::cout << "\nTriangles\n";
    const auto q_tris = AdjacentTriangles<IntT>(e_0, e_1, edge_faces);
    for (const auto t : q_tris) {
        std::cout << t << "  ";
    }
    std::cout << '\n';

    std::cout << "\nInterpolate\n";
    const auto [xi, yi] = Interpolate<IntT>(map, edge_faces, tris, px, py, x, y, x, y);
    std::cout << xi << "  " << yi << '\n';

    return 0;
}