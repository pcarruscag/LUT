rm benchmark.o benchmark
g++ -pedantic -Wextra -Wall -Wconversion -Werror -std=c++17 -m64 -O2 -march=native -funroll-loops -ftree-vectorize -fopenmp -ffast-math -fno-finite-math-only -DNDEBUG -c ./benchmark.cpp -o ./benchmark.o
g++ -o benchmark ./benchmark.o -m64 -s -fopenmp

./benchmark