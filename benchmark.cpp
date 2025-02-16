#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
#include <cmath>

#include "containers.hpp"
#include "trap_map.hpp"


int main() {
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
    auto t_build = -omp_get_wtime();

    Matrix2i edge_pts, edge_faces;
    ExtractEdges<IntT>(tris, edge_pts, edge_faces);

    TrapezoidalMap<VectorInt, VectorReal> map;
    BuildTrapezoidalMap<IntT>(edge_pts, x, y, map);

    t_build += omp_get_wtime();
    std::cout << "Built map in " << t_build << "s\n";

    // Manufacture a function.
    auto func = [](const auto x, const auto y) {
        return sin((x - y + 2) * 0.5 * M_PI);
    };

    VectorReal z;
    z.resize(n_pt);
    for (IntT i = 0; i < n_pt; ++i) z[i] = func(x[i], y[i]);

    // Query the map at triangle centroids.
    VectorReal xi, yi, zi;
    xi.resize(n_tri);
    yi.resize(n_tri);
    zi.resize(n_tri);

    for (IntT i = 0; i < n_tri; ++i) {
        xi[i] = (x[tris(i, 0)] + x[tris(i, 1)] + x[tris(i, 2)]) / 3;
        yi[i] = (y[tris(i, 0)] + y[tris(i, 1)] + y[tris(i, 2)]) / 3;
    }

    auto t_query = -omp_get_wtime();
    for (IntT i = 0; i < n_tri; ++i) {
        zi[i] = Interpolate<IntT>(map, edge_faces, tris, xi[i], yi[i], x, y, z)[0];
    }
    t_query += omp_get_wtime();
    std::cout << "Query map in " << t_query << "s\n";

    RealT max_err = 0, avg_err = 0;
    for (IntT i = 0; i < n_tri; ++i) {
        const auto val = func(xi[i], yi[i]);
        const auto err = fabs(val - zi[i]);
        max_err = fmax(max_err, err);
        avg_err += err;
    }
    avg_err /= static_cast<RealT>(n_tri);
    std::cout << "Max/Avg error: " << max_err << "/" << avg_err << "\n";

    return 0;
}