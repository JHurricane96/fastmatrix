#include <cassert>
#include <vector>

namespace fastmatrix {

template <class E, typename T> class expression {
public:
  E const &get_const_derived() const { return static_cast<E const &>(*this); }

  T operator()(std::size_t i, std::size_t j) const {
    return get_const_derived()(i, j);
  }

  std::size_t num_rows() const { return get_const_derived().num_rows(); }

  std::size_t num_cols() const { return get_const_derived().num_cols(); }
};

template <typename T> class matrix : public expression<matrix<T>, T> {
private:
  std::vector<T> container;
  std::size_t n_rows;
  std::size_t n_cols;

public:
  matrix() {}

  matrix(std::size_t n_rows, std::size_t n_cols)
      : container(n_rows * n_cols), n_rows(n_rows), n_cols(n_cols) {}

  matrix(std::size_t n_rows, std::size_t n_cols, T fill)
      : container(n_rows * n_cols, fill), n_rows(n_rows), n_cols(n_cols) {}

  matrix(matrix<T> &&other)
      : container(std::move(other.container)), n_rows(other.num_rows()),
        n_cols(other.num_cols()) {}

  template <typename E>
  matrix(expression<E, T> const &other)
      : n_rows(other.num_rows()), n_cols(other.num_cols()),
        container(other.num_rows() * other.num_cols()) {
    assign(other.get_const_derived());
  }

  template <typename E> matrix &operator=(expression<E, T> const &other) {
    assert(n_rows >= other.num_rows());
    assert(n_cols >= other.num_cols());
    assign(other.get_const_derived());
    return *this;
  }

  matrix &operator=(matrix<T> &&other) {
    if (this != &other) {
      container = std::move(other.container);
      n_rows = other.n_rows;
      n_cols = other.n_cols;
    }
    return *this;
  }

  T operator()(std::size_t i, std::size_t j) const {
    assert(i < n_rows && i >= 0);
    assert(j < n_cols && j >= 0);
    return container[i * n_rows + j];
  }

  std::vector<T> &get_container() { return container; }

  std::size_t num_rows() const { return n_rows; }

  std::size_t num_cols() const { return n_cols; }

  const matrix<T> &eval() const { return *this; }

  void set_elt(std::size_t i, std::size_t j, T value) {
    container[i * n_rows + j] = value;
  }

  template <typename E, typename ET>
  void assign(expression<E, ET> const &expr) {
    for (std::size_t i = 0; i < n_rows; ++i) {
      for (std::size_t j = 0; j < n_cols; ++j) {
        container[i * n_rows + j] = expr(i, j);
      }
    }
  }
};

template <class E1, class E2, typename T, class Op>
class cwise_matrix_binary_operation
    : public expression<cwise_matrix_binary_operation<E1, E2, T, Op>, T> {
private:
  E1 const &expr1;
  E2 const &expr2;

public:
  cwise_matrix_binary_operation(expression<E1, T> const &expr1,
                                expression<E2, T> const &expr2)
      : expr1(expr1.get_const_derived()), expr2(expr2.get_const_derived()) {}

  T operator()(std::size_t i, std::size_t j) const {
    return Op::apply(expr1, expr2, i, j);
  }

  const matrix<T> eval() const {
    matrix<T> temp(num_rows(), num_cols());
    temp.assign((*this));
    return temp;
  }

  std::size_t num_rows() const { return expr1.num_rows(); }

  std::size_t num_cols() const { return expr1.num_cols(); }
};

template <typename E1, typename E2, typename T>
class matrix_product : public expression<matrix_product<E1, E2, T>, T> {
private:
  E1 const &expr1;
  E2 const &expr2;
  matrix<T> temp;

public:
  matrix_product(expression<E1, T> const &expr1, expression<E2, T> const &expr2)
      : expr1(expr1.get_const_derived()), expr2(expr2.get_const_derived()),
        temp(expr1.num_rows(), expr2.num_cols()) {
    for (std::size_t i = 0; i < num_rows(); ++i) {
      for (std::size_t j = 0; j < num_cols(); ++j) {
        temp.set_elt(i, j, expr1(i, 0) * expr2(0, j));
        for (std::size_t k = 1; k < expr1.num_cols(); ++k) {
          temp.set_elt(i, j, expr1(i, k) * expr2(k, j) + temp(i, j));
        }
      }
    }
  }

  T operator()(std::size_t i, std::size_t j) const { return temp(i, j); }

  matrix<T> const &eval() const { return temp; }

  std::size_t num_rows() const { return expr1.num_rows(); }

  std::size_t num_cols() const { return expr2.num_cols(); }
};

template <class E1, class E2, typename T> struct matrix_add {
  static T apply(expression<E1, T> const &expr1, expression<E2, T> const &expr2,
                 std::size_t i, std::size_t j) {
    return expr1(i, j) + expr2(i, j);
  }
};

template <class E1, class E2, typename T>
cwise_matrix_binary_operation<E1, E2, T, matrix_add<E1, E2, T>>
operator+(expression<E1, T> const &expr1, expression<E2, T> const &expr2) {
  assert(expr1.num_rows() == expr2.num_rows());
  assert(expr1.num_cols() == expr2.num_cols());
  return cwise_matrix_binary_operation<E1, E2, T, matrix_add<E1, E2, T>>(expr1,
                                                                         expr2);
}

template <class E1, class E2, typename T>
matrix_product<E1, E2, T> operator*(expression<E1, T> const &expr1,
                                    expression<E2, T> const &expr2) {
  assert(expr1.num_cols() == expr2.num_rows());
  return matrix_product<E1, E2, T>(expr1, expr2);
}
}
