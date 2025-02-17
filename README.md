Trapezoidal map implementation using STL algorithms.
Sequential and parallel (PSTL + OpenMP) versions of the build algorithm.
It takes ~2 threads for the parallel version to break even with the sequential
implementation (with more threads it becomes much faster).
Querying is the same for both versions, and thread-safe.