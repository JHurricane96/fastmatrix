#include "fastmatrix/fastmatrix.hpp"
#include <iostream>
#include <typeinfo>

using namespace fastmatrix;

int main() {
	matrix<int> a(3, 3, 4);
	matrix<int> b(3, 3, 5);
	matrix<int> c(3, 3, 5);
	matrix<int> d(10, 10, 5);
	auto e = d;
	// std::cout << typeid(a + b + c + d + e).name() << std::endl;
	matrix<int> x = d;
	x = x + x;
	std::cout << x(2, 2) << std::endl;
	return 0;
}