#include <cassert>
#include <vector>

namespace fastmatrix {

template <class E, typename T> class matrix_expression {
public:
  T operator()(std::size_t i, std::size_t j) const {
    return static_cast<E const &>(*this)(i, j);
  }

  std::size_t num_rows() const {
    return static_cast<E const &>(*this).num_rows();
  }

  std::size_t num_cols() const {
    return static_cast<E const &>(*this).num_cols();
  }

  operator E &() { return static_cast<E &>(*this); }
  operator const E &() const { return static_cast<const E &>(*this); }
};

template <typename T> class matrix : public matrix_expression<matrix<T>, T> {
private:
  std::vector<T> container;
  std::size_t n_rows;
  std::size_t n_cols;

public:
  matrix(std::size_t n_rows, std::size_t n_cols)
      : container(n_rows * n_cols), n_rows(n_rows), n_cols(n_cols) {}

  matrix(std::size_t n_rows, std::size_t n_cols, T fill)
      : container(n_rows * n_cols, fill), n_rows(n_rows), n_cols(n_cols) {}

  template <typename E>
  matrix(matrix_expression<E, T> const &other)
      : n_rows(other.num_rows()), n_cols(other.num_cols()),
        container(other.num_rows() * other.num_cols()) {
    for (std::size_t i = 0; i < n_rows; ++i) {
      for (std::size_t j = 0; j < n_cols; ++j) {
        container[i * n_rows + j] = other(i, j);
      }
    }
  }

  T operator()(std::size_t i, std::size_t j) const {
    assert(i < n_rows && i >= 0);
    assert(j < n_cols && j >= 0);
    return container[i * n_rows + j];
  }

  std::vector<T> &get_container() { return container; }

  std::size_t num_rows() const { return n_rows; }

  std::size_t num_cols() const { return n_cols; }
};

template <class E1, class E2, typename T, class Op>
class matrix_binary_operation
    : public matrix_expression<matrix_binary_operation<E1, E2, T, Op>, T> {
public:
  E1 const &expr1;
  E2 const &expr2;

public:
  matrix_binary_operation(E1 const &expr1, E2 const &expr2)
      : expr1(expr1), expr2(expr2) {}

  T operator()(std::size_t i, std::size_t j) const {
    return Op::apply(expr1, expr2, i, j);
  }

  std::size_t num_rows() const { return expr1.num_rows(); }

  std::size_t num_cols() const { return expr1.num_cols(); }
};

template <class E1, class E2, typename T> struct matrix_add {
  static T apply(matrix_expression<E1, T> const &expr1,
                 matrix_expression<E2, T> const &expr2, std::size_t i,
                 std::size_t j) {
    return expr1(i, j) + expr2(i, j);
  }
};

template <class E1, class E2, typename T>
matrix_binary_operation<E1, E2, T, matrix_add<E1, E2, T>>
operator+(matrix_expression<E1, T> const &expr1,
          matrix_expression<E2, T> const &expr2) {
  assert(expr1.num_rows() == expr2.num_rows());
  assert(expr1.num_cols() == expr2.num_cols());
  return matrix_binary_operation<E1, E2, T, matrix_add<E1, E2, T>>(expr1,
                                                                   expr2);
}

template <class E1, class E2, typename T> struct matrix_multiply {
  static T apply(matrix_expression<E1, T> const &expr1,
                 matrix_expression<E2, T> const &expr2, std::size_t i,
                 std::size_t j) {
    T ans = expr1(i, 0) * expr2(0, j);
    for (int k = 1; k < expr1.num_cols(); ++k) {
      ans = ans + expr1(i, k) * expr2(k, j);
    }
    return ans;
  }
};

template <class E1, class E2, typename T>
matrix_binary_operation<E1, E2, T, matrix_multiply<E1, E2, T>>
operator*(matrix_expression<E1, T> const &expr1,
          matrix_expression<E2, T> const &expr2) {
  assert(expr1.num_cols() == expr2.num_rows());
  return matrix_binary_operation<E1, E2, T, matrix_multiply<E1, E2, T>>(expr1,
                                                                        expr2);
}
}