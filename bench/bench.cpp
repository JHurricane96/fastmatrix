#include "fastmatrix/fastmatrix.hpp"
#include "timer.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <random>

using namespace std;
using namespace fastmatrix;

const int TRIES = 10;
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

double square(const double num) {
  return num * num;
}

double timer_mean(const BenchTimer &timer) {
  return timer.total() / TRIES;
}

double timer_sd(const BenchTimer &timer) {
  double mean_squares = timer.squared_total() / TRIES;
  double squared_mean = square(timer_mean(timer));
  return sqrt(mean_squares - squared_mean);
}

void print_results(string test_name, const BenchTimer &timer, double flops_factor) {
  auto mean = timer_mean(timer);
  auto sd = timer_sd(timer);
  auto flops = flops_factor * REPEAT * TRIES / (pow(1024., 3) * timer.total());
  std::cout << test_name << " :: Mean: " << mean << "s, SD: " << sd << "s; " << flops << " GFlops\n";
}

int main() {
  size_t row = 1E4, col = 1E3;
  auto a = make_rand_matrix<float>(row, col);
  auto b = make_rand_matrix<float>(row, col);
  auto c = make_rand_matrix<float>(row, col);
  auto d = make_rand_matrix<float>(row, col);
  auto result = matrix<float>(row, col, 0);

  BenchTimer timer;

  BENCH(timer, TRIES, REPEAT, result = (((a + b).eval() + c).eval() + d));
  print_results("a+b+c+d (Eager)", timer, double(row * col * 3));

  BENCH(timer, TRIES, REPEAT, result = a + b + c + d);
  print_results("a+b+c+d (Lazy)", timer, double(row * col * 3));

  BENCH(timer, TRIES, REPEAT, result = ((a * 5).eval() + (b * 5).eval()));
  print_results("(5*a)+(5*b) (Eager)", timer, double(row * col * 3));

  BENCH(timer, TRIES, REPEAT, result = a * 5 + b * 5);
  print_results("(5*a)+(5*b) (Lazy)", timer, double(row * col * 3));

  return 0;
}
