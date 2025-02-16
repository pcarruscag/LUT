#include <cstddef>
#include <cstdint>
#include <vector>

template <typename T, size_t N>
struct Matrix {
    std::vector<T> data;

    void resize(size_t rows, size_t) { data.resize(rows * N); }
    const auto rows() const { return data.size() / N; }

    const auto& operator()(size_t i, size_t j) const { return data[i * N + j]; }
    auto& operator()(size_t i, size_t j) { return data[i * N + j]; }
};

using IntT = int32_t;
using RealT = double;

using Matrix2i = Matrix<IntT, 2>;
using Matrix3i = Matrix<IntT, 3>;

using VectorInt = std::vector<IntT>;
using VectorReal = std::vector<RealT>;