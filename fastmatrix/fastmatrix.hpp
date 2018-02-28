#include <cassert>
#include <ostream>
#include <type_traits>
#include <vector>

namespace fastmatrix {

template <typename T>
struct eval_return_type {
  using type = typename T::EvalReturnType;
};

template <typename T>
using eval_return_type_t = typename eval_return_type<T>::type;

template <typename T>
struct element_type {
  using type = typename T::ElementType;
};

template <typename T>
using element_type_t = typename element_type<T>::type;

template <typename T>
struct storage_type {
  using type = T const &;
};

template <typename T>
using storage_type_t = typename storage_type<T>::type;

template <typename E>
class expression {
public:
  inline E const &get_const_derived() const {
    return static_cast<E const &>(*this);
  }

  inline std::size_t num_rows() const {
    return get_const_derived().num_rows();
  }

  inline std::size_t num_cols() const {
    return get_const_derived().num_cols();
  }
};

template <typename T>
class scalar_expression : public expression<scalar_expression<T>> {
private:
  T scalar;

public:
  using EvalReturnType = T;
  using ElementType = T;

  inline scalar_expression(T scalar) : scalar(scalar) {}

  inline T operator()(std::size_t, std::size_t) const {
    return scalar;
  }

  inline T eval() const {
    return scalar;
  }

  inline std::size_t num_rows() const {
    return 0;
  }

  inline std::size_t num_cols() const {
    return 0;
  }
};

template <typename T>
class matrix : public expression<matrix<T>> {
private:
  std::vector<T> container;
  std::size_t n_rows;
  std::size_t n_cols;

public:
  using EvalReturnType = matrix<T>;
  using ElementType = T;

  inline matrix() {}

  inline matrix(std::size_t n_rows, std::size_t n_cols)
      : container(n_rows * n_cols), n_rows(n_rows), n_cols(n_cols) {}

  inline matrix(std::size_t n_rows, std::size_t n_cols, T fill)
      : container(n_rows * n_cols, fill), n_rows(n_rows), n_cols(n_cols) {}

  inline matrix(matrix<T> &&other)
      : container(std::move(other.container)), n_rows(other.num_rows()), n_cols(other.num_cols()) {}

  template <typename E>
  inline matrix(expression<E> const &other)
      : container(other.num_rows() * other.num_cols()), n_rows(other.num_rows()),
        n_cols(other.num_cols()) {
    assign(other.get_const_derived());
  }

  template <typename E>
  inline matrix &operator=(expression<E> const &other) {
    assert(n_rows >= other.num_rows());
    assert(n_cols >= other.num_cols());
    assign(other.get_const_derived());
    return *this;
  }

  inline matrix &operator=(matrix<T> &&other) {
    if (this != &other) {
      container = std::move(other.container);
      n_rows = other.n_rows;
      n_cols = other.n_cols;
    }
    return *this;
  }

  inline T operator()(std::size_t i, std::size_t j) const {
    assert(i < n_rows && i >= 0);
    assert(j < n_cols && j >= 0);
    return container[i * n_cols + j];
  }

  inline std::vector<T> &get_container() {
    return container;
  }

  inline std::size_t num_rows() const {
    return n_rows;
  }

  inline std::size_t num_cols() const {
    return n_cols;
  }

  inline const matrix<T> &eval() const {
    return *this;
  }

  inline void set_elt(std::size_t i, std::size_t j, T value) {
    container[i * n_cols + j] = value;
  }

  template <typename E>
  inline void assign(expression<E> const &expr) {
    for (std::size_t i = 0; i < n_rows; ++i) {
      for (std::size_t j = 0; j < n_cols; ++j) {
        container[i * n_cols + j] = expr.get_const_derived()(i, j);
      }
    }
  }

  friend std::ostream &operator<<(std::ostream &ostream, const matrix<T> &mat) {
    for (std::size_t i = 0; i < mat.n_rows; ++i) {
      for (std::size_t j = 0; j < mat.n_cols; ++j) {
        ostream << mat(i, j) << ", ";
      }
      ostream << "\n";
    }
    return ostream;
  }

  template <typename E>
  inline matrix<T> &operator+=(expression<E> const &expr);

  template <typename E>
  inline matrix<T> &operator*=(expression<E> const &expr);

  template <typename E>
  inline matrix<T> &operator-=(expression<E> const &expr);

  template <typename Scalar, typename = std::enable_if_t<std::is_arithmetic<Scalar>::value>>
  inline matrix<T> &operator+=(Scalar const &expr);

  template <typename Scalar, typename = std::enable_if_t<std::is_arithmetic<Scalar>::value>>
  inline matrix<T> &operator-=(Scalar const &expr);

  template <typename Scalar, typename = std::enable_if_t<std::is_arithmetic<Scalar>::value>>
  inline matrix<T> &operator*=(Scalar const &expr);
};
} // namespace fastmatrix

namespace std {
using namespace fastmatrix;

template <typename T1, typename T2>
struct common_type<matrix<T1>, matrix<T2>> {
  using type = matrix<std::common_type_t<T1, T2>>;
};

template <typename T1, typename T2>
struct common_type<matrix<T1>, T2> {
  using type = matrix<std::common_type_t<T1, T2>>;
};
} // namespace std

namespace fastmatrix {

template <typename T>
struct storage_type<scalar_expression<T>> {
  using type = scalar_expression<T>;
};

template <typename Op, typename E1, typename E2>
class cwise_matrix_binary_operation : public expression<cwise_matrix_binary_operation<Op, E1, E2>> {
private:
  E1 const &expr1;
  storage_type_t<E2> expr2;

public:
  using EvalReturnType = std::common_type_t<eval_return_type_t<E1>, eval_return_type_t<E2>>;
  using ElementType = std::common_type_t<element_type_t<E1>, element_type_t<E2>>;

  inline cwise_matrix_binary_operation(expression<E1> const &expr1, expression<E2> const &expr2)
      : expr1(expr1.get_const_derived()), expr2(expr2.get_const_derived()) {}

  inline ElementType operator()(std::size_t i, std::size_t j) const {
    return Op::apply(expr1, expr2, i, j);
  }

  inline const EvalReturnType eval() const {
    EvalReturnType temp(num_rows(), num_cols());
    temp.assign((*this));
    return temp;
  }

  inline std::size_t num_rows() const {
    return expr1.num_rows();
  }

  inline std::size_t num_cols() const {
    return expr1.num_cols();
  }
};

template <template <typename E1, typename E2> typename Op, typename E1, typename E2>
cwise_matrix_binary_operation<Op<E1, E2>, E1, E2>
make_cwise_matrix_binary_operation(expression<E1> const &expr1, expression<E2> const &expr2) {
  return cwise_matrix_binary_operation<Op<E1, E2>, E1, E2>(expr1, expr2);
}

template <typename E1, typename E2>
class matrix_product : public expression<matrix_product<E1, E2>> {
public:
  using EvalReturnType = std::common_type_t<eval_return_type_t<E1>, eval_return_type_t<E2>>;
  using ElementType = std::common_type_t<element_type_t<E1>, element_type_t<E2>>;

private:
  E1 const &expr1;
  E2 const &expr2;
  EvalReturnType temp;

public:
  inline matrix_product(expression<E1> const &expr1, expression<E2> const &expr2)
      : expr1(expr1.get_const_derived()), expr2(expr2.get_const_derived()),
        temp(expr1.num_rows(), expr2.num_cols()) {
    for (std::size_t i = 0; i < num_rows(); ++i) {
      for (std::size_t j = 0; j < num_cols(); ++j) {
        temp.set_elt(i, j, expr1.get_const_derived()(i, 0) * expr2.get_const_derived()(0, j));
        for (std::size_t k = 1; k < expr1.num_cols(); ++k) {
          temp.set_elt(
              i, j, expr1.get_const_derived()(i, k) * expr2.get_const_derived()(k, j) + temp(i, j));
        }
      }
    }
  }

  inline ElementType operator()(std::size_t i, std::size_t j) const {
    return temp(i, j);
  }

  inline EvalReturnType const &eval() const {
    return temp;
  }

  inline std::size_t num_rows() const {
    return expr1.num_rows();
  }

  inline std::size_t num_cols() const {
    return expr2.num_cols();
  }
};

template <typename E1, typename E2>
struct cwise_matrix_add {
  inline static typename std::common_type_t<element_type_t<E1>, element_type_t<E2>>
  apply(expression<E1> const &expr1, expression<E2> const &expr2, std::size_t i, std::size_t j) {
    return expr1.get_const_derived()(i, j) + expr2.get_const_derived()(i, j);
  }
};

template <typename E1, typename E2>
struct cwise_matrix_multiply {
  inline static typename std::common_type_t<element_type_t<E1>, element_type_t<E2>>
  apply(expression<E1> const &expr1, expression<E2> const &expr2, std::size_t i, std::size_t j) {
    return expr1.get_const_derived()(i, j) * expr2.get_const_derived()(i, j);
  }
};

template <typename E1, typename E2>
struct cwise_matrix_subtract {
  inline static typename std::common_type_t<element_type_t<E1>, element_type_t<E2>>
  apply(expression<E1> const &expr1, expression<E2> const &expr2, std::size_t i, std::size_t j) {
    return expr1.get_const_derived()(i, j) - expr2.get_const_derived()(i, j);
  }
};

template <typename E1, typename E2>
inline auto operator+(expression<E1> const &expr1, expression<E2> const &expr2) {
  assert(expr1.num_rows() == expr2.num_rows());
  assert(expr1.num_cols() == expr2.num_cols());
  return make_cwise_matrix_binary_operation<cwise_matrix_add>(expr1, expr2);
}

template <typename E, typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value, T>>
inline auto operator+(expression<E> const &expr, T const &scalar) {
  return make_cwise_matrix_binary_operation<cwise_matrix_add>(expr, scalar_expression(scalar));
}

template <typename E1, typename E2>
inline auto operator*(expression<E1> const &expr1, expression<E2> const &expr2) {
  assert(expr1.num_cols() == expr2.num_rows());
  return matrix_product(expr1, expr2);
}

template <typename E, typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value, T>>
inline auto operator*(expression<E> const &expr, T const &scalar) {
  return make_cwise_matrix_binary_operation<cwise_matrix_multiply>(expr, scalar_expression(scalar));
}

template <typename E1, typename E2>
inline auto operator-(expression<E1> const &expr1, expression<E2> const &expr2) {
  assert(expr1.num_rows() == expr2.num_rows());
  assert(expr1.num_cols() == expr2.num_cols());
  return make_cwise_matrix_binary_operation<cwise_matrix_subtract>(expr1, expr2);
}

template <typename E, typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value, T>>
inline auto operator-(expression<E> const &expr, T const &scalar) {
  return make_cwise_matrix_binary_operation<cwise_matrix_subtract>(expr, scalar_expression(scalar));
}

template <typename T>
template <typename E>
inline matrix<T> &matrix<T>::operator+=(expression<E> const &expr) {
  assert(n_rows == expr.num_rows());
  assert(n_cols == expr.num_cols());
  assign(make_cwise_matrix_binary_operation<cwise_matrix_add>(*this, expr));
  return *this;
}

template <typename T>
template <typename Scalar, typename>
inline matrix<T> &matrix<T>::operator+=(Scalar const &scalar) {
  assign(make_cwise_matrix_binary_operation<cwise_matrix_add>(*this, scalar_expression(scalar)));
  return *this;
}

template <typename T>
template <typename E>
inline matrix<T> &matrix<T>::operator*=(expression<E> const &expr) {
  assert(n_rows == expr.num_cols());
  assign(matrix_product(*this, expr));
  return *this;
}

template <typename T>
template <typename Scalar, typename>
inline matrix<T> &matrix<T>::operator*=(Scalar const &scalar) {
  assign(
      make_cwise_matrix_binary_operation<cwise_matrix_multiply>(*this, scalar_expression(scalar)));
  return *this;
}

template <typename T>
template <typename E>
inline matrix<T> &matrix<T>::operator-=(expression<E> const &expr) {
  assert(n_rows == expr.num_rows());
  assert(n_cols == expr.num_cols());
  assign(make_cwise_matrix_binary_operation<cwise_matrix_subtract>(*this, expr));
  return *this;
}

template <typename T>
template <typename Scalar, typename>
inline matrix<T> &matrix<T>::operator-=(Scalar const &scalar) {
  assign(
      make_cwise_matrix_binary_operation<cwise_matrix_subtract>(*this, scalar_expression(scalar)));
  return *this;
}
} // namespace fastmatrix
