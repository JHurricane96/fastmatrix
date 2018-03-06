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
matrix<T> make_rand_matrix(size_t num_rows, size_t num_cols) {
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
  size_t row = 1E2, col = 1E2;
  auto a = make_rand_matrix<float>(row, col);
  auto b = make_rand_matrix<float>(row, col);
  auto c = make_rand_matrix<float>(row, col);
  auto d = make_rand_matrix<float>(row, col);
  auto result = matrix<float>(row, col, 0);

  BenchTimer timer;

  BENCH(timer, TRIES, REPEAT, result = (((a + b).eval() + c).eval() + d));
  print_results("a+b+c+d (Eager)", timer.value(), (double(row * col * REPEAT * 3) / timer.value()));

  BENCH(timer, TRIES, REPEAT, result = a + b + c + d);
  print_results("a+b+c+d (Lazy)", timer.value(), (double(row * col * REPEAT * 3) / timer.value()));

  BENCH(timer, TRIES, REPEAT, result = a * b * c * d);
  print_results("a*b*c*d", timer.value(), (double(row * col * col * REPEAT * 3) / timer.value()));

  BENCH(timer, TRIES, REPEAT, result = a * 5 + b);
  print_results("(5*a)+b", timer.value(), (double(row * col * REPEAT * 2) / timer.value()));

  BENCH(timer, TRIES, REPEAT, result = result + b * 5 - 10 - a * result;);
  print_results("result + b*5 - 10 - a*result", timer.value(),
                (double(row * col * (col + 4) * REPEAT) / timer.value()));
  return 0;
}
