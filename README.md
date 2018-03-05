# fastmatrix

A linear algebra library in C++ to do fast math.

Uses [expression templates](https://en.wikipedia.org/wiki/Expression_templates) to build expression trees at compile time. Most expressions are lazily evaluated.

## Usage

fastmatrix is header-only, simply add the location of fastmatrix.hpp to your include path while compiling.

For optimal speed, compile with the highest level of optimizations (`-O3` on gcc/clang) and with assertions disabled (`-DNDEBUG` on gcc/clang)

Example code:

```cpp
#include "fastmatrix/fastmatrix.hpp"
#include <iostream>

using namespace fastmatrix;

int main () {
	// Make a 3x3 matrix filled with 0.5
	matrix<float> a(3, 3, 0.5); // (row, col, value)
	matrix<float> b(3, 3, 0.4);

	// Make an empty 3x3 matrix
	matrix<float> result(3, 3); // (row, col)

	// Basic arithmetic can be done easily
	result = a + b;
	result = a - b;
	result = a * b;

	// Element-wise arithmetic can be done too
	// Here, 5 is subtracted from every element of a
	result = a - 5;

	// Expressions can be as long as you want them to be
	result = 5 * a - result + b * 10.5 - result;

	// Use eval to trigger eager evaluation.
	// Here, it is beneficial to calculate (b + result)
	// in advance to prevent repeated calculations
	result = a * (b + result).eval();

	// Getting and setting elements is possible
	// All indexing is 0-based
	std::cout << result(0, 1) << '\n'; // (row, column)
	result.set_elt(0, 1, -1.5); // (row, column, new value)
	std::cout << result(0, 1) << '\n';

	// Printing a matrix is easy
	std::cout << result << '\n';
	return 0;
}
```

## Development

To run the tests and benchmarks on Linux:

(Dependencies: CMake >= 3.10.1, g++ >= 7/clang++ >= 5)

1. Clone this repo

2. `mkdir build && cd build`

3. `cmake ..` (Add `-DBUILD_TESTS=OFF` to not build tests, `-DBUILD_BENCHMARKS=OFF` to not build benchmarks). You can set the compiler to use by setting the `CXX` environment variable before running cmake.

4. `cmake --build .`

5. `./tests` to run the tests and `./bench/bench` to run the benchmarks

Run `doxygen` to generate the documentation in the `/docs` folder.
