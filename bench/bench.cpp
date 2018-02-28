#include "fastmatrix/fastmatrix.hpp"
#include "timer.h"
#include <iostream>
#include <limits>
#include <random>

using namespace std;
using namespace fastmatrix;

const int TRIES = 1;
const int REPEAT = 100;

template <typename T>
matrix<T> make_matrix(size_t num_rows, size_t num_cols) {
  static_assert(is_floating_point<T>::value);
  static std::mt19937 gen{std::random_device{}()};
  static std::uniform_real_distribution dis(numeric_limits<T>::min(), numeric_limits<T>::max());

  matrix<T> m(num_rows, num_cols);
  for (size_t i = 0; i < num_rows; ++i) {
    for (size_t j = 0; j < num_cols; ++j) {
      m.set_elt(i, j, dis(gen));
    }
  }
  return m;
}

void print_results(string test_name, double time_s, double flops) {
  std::cout << test_name << ": " << time_s << "s, " << flops / (1024. * 1024. * 1024.)
            << " GFlops\n";
}

int main() {
  size_t row = 1E3, col = 1E3;
  BenchTimer timer;
  timer.reset();
  timer.start();
  auto a = make_matrix<float>(row, col);
  auto b = make_matrix<float>(row, col);
  auto c = make_matrix<float>(row, col);
  auto d = make_matrix<float>(row, col);
  auto result = matrix<float>(row, col, 0);
  timer.stop();
  std::cout << timer.value() << std::endl;
  BENCH(timer, TRIES, REPEAT, result = (((a + b).eval() + c).eval() + d));
  print_results("Eager Add", timer.value(), (double(row * col * REPEAT * 3) / timer.value()));
  BENCH(timer, TRIES, REPEAT, result = a + b + c + d);
  print_results("Lazy Add", timer.value(), (double(row * col * REPEAT * 3) / timer.value()));
  return 0;
}