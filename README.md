Trapezoidal map implementation using STL algorithms.
Sequential and parallel (PSTL + OpenMP) versions of the build algorithm.
With a single thread, the parallel version is ~25% slower than the sequential
implementation (with more threads it becomes much faster).
Querying is the same for both versions, and thread-safe.