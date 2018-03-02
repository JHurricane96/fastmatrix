# fastmatrix

A linear algebra library in C++ to do fast math.

## Usage

fastmatrix is header-only, simply add the location of fastmatrix.hpp to your include path while compiling.

## Development

To run the tests and benchmarks on Linux:

(Dependencies: CMake >= 3.10.1, g++ >= 7/clang++ >= 5)

1. Clone this repo

2. `mkdir build && cd build`

3. `cmake ..` (Add -DBUILD_TESTS=OFF to not build tests, -DBUILD_BENCHMARKS=OFF to not build benchmarks)

4. `make`

5. `./tests` to run the tests and `./bench/bench` to run the benchmarks

Run `doxygen` to generate the documentation in the `/docs` folder.
