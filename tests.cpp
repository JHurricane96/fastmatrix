#include "fastmatrix/fastmatrix.hpp"
#include <iostream>

using namespace std;

template <typename T, typename F>
matrix<T> make_matrix(size_t num_rows, size_t num_cols, F const &filler) {
  matrix<T> m(num_rows, num_cols);
  for (size_t i = 0; i < num_rows; ++i) {
    for (size_t j = 0; j < num_cols; ++j) {
      m.set_elt(i, j, filler(i, j));
    }
  }
  return m;
}

template <typename T, typename F>
vector<vector<T>> make_vector(size_t num_rows, size_t num_cols, F const &filler) {
  vector<vector<T>> v(num_rows);
  for (size_t i = 0; i < v.size(); ++i) {
    auto &row = v[i];
    for (size_t j = 0; j < num_cols; ++j) {
      row.push_back(filler(i, j));
    }
  }
  return v;
}

template <typename T>
void print_vector(vector<vector<T>> const &v) {
  for (auto &row : v) {
    for (auto &elt : row) {
      std::cout << elt << ", ";
    }
    std::cout << '\n';
  }
}

template <typename T>
void assert_equals(vector<vector<T>> const &v, matrix<T> const &m) {
  for (size_t i = 0; i < m.num_rows(); ++i) {
    for (size_t j = 0; j < m.num_cols(); ++j) {
      assert(v[i][j] == m(i, j));
    }
  }
}

void test_matrix_add() {
  size_t row = 1000, col = 3000;

  auto a_filler = [col](size_t i, size_t j) { return i * col + j; };
  auto b_filler = [row, col](size_t i, size_t j) { return i * col + j + row * col; };
  auto result_filler = [&a_filler, &b_filler](size_t i, size_t j) {
    return a_filler(i, j) + b_filler(i, j);
  };
  auto empty_filler = [](size_t, size_t) { return 0; };

  auto v_a = make_vector<int>(row, col, a_filler);
  auto v_b = make_vector<int>(row, col, b_filler);
  auto v_result = make_vector<int>(row, col, result_filler);

  auto m_a = make_matrix<int>(row, col, a_filler);
  auto m_b = make_matrix<int>(row, col, b_filler);
  auto m_result = make_matrix<int>(row, col, empty_filler);

  m_result = m_a + m_b;

  assert_equals(v_result, m_result);
  std::cout << "Successfully executed test_matrix_add\n";
}

void test_matrix_multiply() {
  size_t row = 100, col = 300, mid_row = 1200;

  auto a_filler = [mid_row](size_t i, size_t j) { return i * mid_row + j; };
  auto b_filler = [mid_row, col](size_t i, size_t j) { return i * col + j + mid_row * col; };
  auto result_filler = [&a_filler, &b_filler, mid_row](size_t i, size_t j) {
    size_t result = 0;
    for (size_t k = 0; k < mid_row; ++k) {
      result += a_filler(i, k) * b_filler(k, j);
    }
    return result;
  };
  auto empty_filler = [](size_t, size_t) { return 0; };

  auto v_a = make_vector<int>(row, mid_row, a_filler);
  auto v_b = make_vector<int>(mid_row, col, b_filler);
  auto v_result = make_vector<int>(row, col, result_filler);

  auto m_a = make_matrix<int>(row, mid_row, a_filler);
  auto m_b = make_matrix<int>(mid_row, col, b_filler);
  auto m_result = make_matrix<int>(row, col, empty_filler);

  m_result = m_a * m_b;

  assert_equals(v_result, m_result);
  std::cout << "Successfully executed test_matrix_multiply\n";
}

void test_matrix_subtract() {
  size_t row = 1000, col = 3000;

  auto a_filler = [col](size_t i, size_t j) { return i * col + j; };
  auto b_filler = [row, col](size_t i, size_t j) { return i * col + j + row * col; };
  auto result_filler = [&a_filler, &b_filler](size_t i, size_t j) {
    return a_filler(i, j) - b_filler(i, j);
  };
  auto empty_filler = [](size_t, size_t) { return 0; };

  auto v_a = make_vector<int>(row, col, a_filler);
  auto v_b = make_vector<int>(row, col, b_filler);
  auto v_result = make_vector<int>(row, col, result_filler);

  auto m_a = make_matrix<int>(row, col, a_filler);
  auto m_b = make_matrix<int>(row, col, b_filler);
  auto m_result = make_matrix<int>(row, col, empty_filler);

  m_result = m_a - m_b;

  assert_equals(v_result, m_result);
  std::cout << "Successfully executed test_matrix_subtract\n";
}

void test_scalar_add() {
  size_t row = 1000, col = 3000;
  const int scalar = 40198;

  auto filler = [col](size_t i, size_t j) { return i * col + j; };
  auto result_filler = [&filler](size_t i, size_t j) { return filler(i, j) + scalar; };
  auto empty_filler = [](size_t, size_t) { return 0; };

  auto v = make_vector<int>(row, col, filler);
  auto v_result = make_vector<int>(row, col, result_filler);

  auto m = make_matrix<int>(row, col, filler);
  auto m_result = make_matrix<int>(row, col, empty_filler);

  m_result = m + scalar;

  assert_equals(v_result, m_result);
  std::cout << "Successfully executed test_scalar_add\n";
}

void test_scalar_multiply() {
  size_t row = 1000, col = 3000;
  const int scalar = 2938;

  auto filler = [col](size_t i, size_t j) { return i * col + j; };
  auto result_filler = [&filler](size_t i, size_t j) { return filler(i, j) * scalar; };
  auto empty_filler = [](size_t, size_t) { return 0; };

  auto v = make_vector<int>(row, col, filler);
  auto v_result = make_vector<int>(row, col, result_filler);

  auto m = make_matrix<int>(row, col, filler);
  auto m_result = make_matrix<int>(row, col, empty_filler);

  m_result = m * scalar;

  assert_equals(v_result, m_result);
  std::cout << "Successfully executed test_scalar_multiply\n";
}

void test_scalar_subtract() {
  size_t row = 1000, col = 3000;
  const int scalar = 2938;

  auto filler = [col](size_t i, size_t j) { return i * col + j; };
  auto result_filler = [&filler](size_t i, size_t j) { return filler(i, j) - scalar; };
  auto empty_filler = [](size_t, size_t) { return 0; };

  auto v = make_vector<int>(row, col, filler);
  auto v_result = make_vector<int>(row, col, result_filler);

  auto m = make_matrix<int>(row, col, filler);
  auto m_result = make_matrix<int>(row, col, empty_filler);

  m_result = m - scalar;

  assert_equals(v_result, m_result);
  std::cout << "Successfully executed test_scalar_subtract\n";
}

void test_variety() {
  size_t row = 5, col = 3;

  auto filler = [col](size_t i, size_t j) { return i * col + j; };
  auto a_filler = [row](size_t i, size_t j) { return i * row + j; };

  int scalar_1 = 5;
  int scalar_2 = 10;
  auto x = make_matrix<int>(row, col, filler);
  auto b = make_matrix<int>(row, col, filler);
  auto a = make_matrix<int>(row, row, a_filler);

  vector<vector<int>> result = {{-100, -104, -108},
                                {-232, -261, -290},
                                {-364, -418, -472},
                                {-496, -575, -654},
                                {-628, -732, -836}};

  x = x + b * scalar_1 - scalar_2 - a * x;
  assert_equals(result, x);
  std::cout << "Successfully executed test_variety\n";
}

int main() {
  std::cout << "Running tests...\n";
  test_matrix_add();
  test_matrix_multiply();
  test_matrix_subtract();
  test_scalar_add();
  test_scalar_multiply();
  test_scalar_subtract();
  test_variety();
  std::cout << "All tests successful!\n";
  return 0;
}