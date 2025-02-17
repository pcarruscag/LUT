rm benchmark.o benchmark
g++ -pedantic -Wextra -Wall -Wconversion -Werror -std=c++17 -m64 -DNDEBUG \
    -O2 -march=native -funroll-loops -ftree-vectorize -ffast-math -fno-finite-math-only \
    -fopenmp -ltbb -c ./benchmark.cpp -o ./benchmark.o
g++ -o benchmark ./benchmark.o -m64 -s -fopenmp -ltbb

./benchmark