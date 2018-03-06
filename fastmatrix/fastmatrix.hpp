#ifndef FASTMATRIX_FASTMATRIX_HPP
#define FASTMATRIX_FASTMATRIX_HPP

#include <cassert>
#include <ostream>
#include <type_traits>
#include <vector>

namespace fastmatrix {

/**
 * \brief      Trait returning the type of the value on calling .eval() on an expression
 *
 * \tparam     T     Expression type
 */
template <typename T>
struct eval_return_type {
  using type = typename T::EvalReturnType;
};

/**
 * Convenience typedef of eval_return_type
 */
template <typename T>
using eval_return_type_t = typename eval_return_type<T>::type;

/**
 * \brief      Trait returning the type of an element in an matrix expression, or the type of a
 * scalar expression
 *
 * \tparam     T     Expression type
 */
template <typename T>
struct element_type {
  using type = typename T::ElementType;
};

/**
 * Convenience typedef of element_type
 */
template <typename T>
using element_type_t = typename element_type<T>::type;

/**
 * \brief      Trait for how to store an expression as another expression's LHS or RHS (whether to
 * store a ref, const ref or by value)
 *
 * \tparam     T     Expression type
 */
template <typename T>
struct storage_type {
  using type = T const &;
};

/**
 * Convenience typedef of storage_type
 */
template <typename T>
using storage_type_t = typename storage_type<T>::type;

/**
 * \brief      Base class for expressions
 *
 * All classes that can be a part of an expression or form an expression themselves inherit from
 * this class. This class uses the Curiously Recurring Template Pattern (CRTP) to enforce static
 * time polymorphism
 *
 * Expressions usually store merely the representation of an expression and not the result of the
 * expression itself to benefit from lazy evaluation and compiler optimizations
 *
 * \tparam     E     Derived class type
 */
template <typename E>
class expression {
public:
  /**
   * \brief      Returns const reference to this object cast to its derived class
   *
   * \return     The cast object
   */
  inline E const &get_const_derived() const {
    return static_cast<E const &>(*this);
  }

  /**
   * \brief      Get number of rows in this matrix expression, 0 if this is a scalar
   *
   * \return     Number of rows
   */
  inline std::size_t num_rows() const {
    return get_const_derived().num_rows();
  }

  /**
   * \brief      Get number of columns in this matrix expression, 0 if this is a scalar
   *
   * \return     Number of columns
   */
  inline std::size_t num_cols() const {
    return get_const_derived().num_cols();
  }
};

/**
 * Trait to determine if a type is not a subclass of expression
 */
template <typename T>
using enable_if_not_expression = std::enable_if_t<!std::is_base_of<expression<T>, T>::value>;

/**
 * \brief      Wrapper expression for scalar values
 *
 * \tparam     T     Type of scalar being wrapped
 */
template <typename T>
class scalar_expression : public expression<scalar_expression<T>> {
private:
  /**
   * Scalar being wrapped by this expression
   */
  T scalar;

public:
  /**
   * Return type of this class's eval method
   */
  using EvalReturnType = T;

  /**
   * Element type of this expression, which is simply the type of the scalar itself
   */
  using ElementType = T;

  /**
   * \brief      Constructor
   *
   * \param[in]  scalar  The scalar to wrap
   */
  inline scalar_expression(T scalar) : scalar(scalar) {}

  /**
   * \brief      Function operator to index into the expression and get an element. Since this is a
   * scalar expression, this simply returns the scalar value
   *
   * \return     The scalar value
   */
  inline T operator()(std::size_t, std::size_t) const {
    return scalar;
  }

  /**
   * \brief      Evaluates this expression
   *
   * Since this is a scalar expression, simply returns the scalar value
   *
   * \return     The scalar value
   */
  inline T eval() const {
    return scalar;
  }

  /**
   * \brief      Number of rows in this expression
   *
   * Since this is a scalar expression, number of rows is 0
   *
   * \return     0
   */
  inline std::size_t num_rows() const {
    return 0;
  }

  /**
   * \brief      Number of columns in this expression
   *
   * Since this is a scalar expression, number of rows is 0
   *
   * \return     0
   */
  inline std::size_t num_cols() const {
    return 0;
  }
};

/**
 * \brief      Class for an actual matrix
 *
 * This is a 2D matrix that stores elements of any scalar type in a container
 *
 * \tparam     T     The type of an element stored in the matrix
 */
template <typename T>
class matrix : public expression<matrix<T>> {
private:
  /**
   * Container in which the matrix elements are stored
   */
  std::vector<T> container;

  /**
   * Number of rows of this matrix
   */
  std::size_t n_rows;

  /**
   * Number of columns of this matrix
   */
  std::size_t n_cols;

public:
  /**
   * Return type of eval() method
   */
  using EvalReturnType = matrix<T>;

  /**
   * Type of elements of this matrix
   */
  using ElementType = T;

  /**
   * \brief      Default constructor to allow empty matrix construction
   */
  inline matrix() = default;

  /**
   * \brief      Constructor
   *
   * Sizes the container appropriately
   *
   * \param[in]  n_rows  The number of rows the matrix should have
   * \param[in]  n_cols  The number of columns the matrix should have
   */
  inline matrix(std::size_t n_rows, std::size_t n_cols)
      : container(n_rows * n_cols), n_rows(n_rows), n_cols(n_cols) {}

  /**
   * \brief      Constructor
   *
   * Fills the container with the given value
   *
   * \param[in]  n_rows  The number of rows the matrix should have
   * \param[in]  n_cols  The number of columns the matrix should have
   * \param[in]  fill    The element with which to fill the container
   */
  inline matrix(std::size_t n_rows, std::size_t n_cols, T fill)
      : container(n_rows * n_cols, fill), n_rows(n_rows), n_cols(n_cols) {}

  /**
   * \brief      Move constructor
   *
   * \param[in]  other  The matrix to construct from
   */
  inline matrix(matrix<T> &&other)
      : container(std::move(other.container)), n_rows(other.num_rows()), n_cols(other.num_cols()) {}

  /**
   * \brief      Constructor from another expression
   *
   * Construct the matrix from the given expression. After building an expression tree, this (and
   * the assignment operator) trigger the actual evaluation of the expression
   *
   * \param      other  The expression
   *
   * \tparam     E      The type of the expression
   */
  template <typename E>
  inline matrix(expression<E> const &other)
      : container(other.num_rows() * other.num_cols()), n_rows(other.num_rows()),
        n_cols(other.num_cols()) {
    assign(other.get_const_derived());
  }

  /**
   * \brief      Assignment from another expression
   *
   * This and the move constructor trigger evaluation of an expression
   *
   * \param      other  The expression
   *
   * \tparam     E      The type of the expression
   *
   * \return     This matrix
   */
  template <typename E>
  inline matrix &operator=(expression<E> const &other) {
    assert(n_rows >= other.num_rows());
    assert(n_cols >= other.num_cols());
    assign(other.get_const_derived());
    return *this;
  }

  /**
   * \brief      Move assignment operator
   *
   * \param[in]  other  The other matrix
   *
   * \return     This matrix
   */
  inline matrix &operator=(matrix<T> &&other) {
    if (this != &other) {
      container = std::move(other.container);
      n_rows = other.n_rows;
      n_cols = other.n_cols;
    }
    return *this;
  }

  /**
   * \brief      Function operator to return elements of this matrix
   *
   * \param[in]  i     Row number of the element to return
   * \param[in]  j     Column number of the element to return
   *
   * \return     The desired element
   */
  inline T operator()(std::size_t i, std::size_t j) const {
    assert(i < n_rows && i >= 0);
    assert(j < n_cols && j >= 0);
    return container[i * n_cols + j];
  }

  /**
   * \brief      Get a reference to the underlying storage container
   *
   * \return     The container
   */
  inline std::vector<T> &get_container() {
    return container;
  }

  /**
   * \brief      Gets number of rows in this matrix
   *
   * \return     Number of rows
   */
  inline std::size_t num_rows() const {
    return n_rows;
  }

  /**
   * \brief      Gets number of columns in this matrix
   *
   * \return     Number of columns
   */
  inline std::size_t num_cols() const {
    return n_cols;
  }

  /**
   * \brief      Evaluate this expression
   *
   * Since this is a matrix whose actual value is already known, simply return a const reference to
   * itself
   *
   * \return     Const reference to this matrix
   */
  inline const matrix<T> &eval() const {
    return *this;
  }

  /**
   * \brief      Set an element of this matrix
   *
   * \param[in]  i      Row number of the element to set
   * \param[in]  j      Column number of this element to set
   * \param[in]  value  The new value of the element
   */
  inline void set_elt(std::size_t i, std::size_t j, T value) {
    container[i * n_cols + j] = value;
  }

  /**
   * \brief      Assign an expression to this matrix
   *
   * \param      expr  The expression to assign
   *
   * \tparam     E     The type of the expression
   */
  template <typename E>
  inline void assign(expression<E> const &expr) {
    for (std::size_t i = 0; i < n_rows; ++i) {
      for (std::size_t j = 0; j < n_cols; ++j) {
        container[i * n_cols + j] = expr.get_const_derived()(i, j);
      }
    }
  }

  /**
   * \brief      The stream operator to print the matrix easily
   *
   * \param      ostream  The output stream
   * \param[in]  mat      The matrix
   *
   * \return     The output stream
   */
  friend std::ostream &operator<<(std::ostream &ostream, const matrix<T> &mat) {
    for (std::size_t i = 0; i < mat.n_rows; ++i) {
      for (std::size_t j = 0; j < mat.n_cols; ++j) {
        ostream << mat(i, j) << ", ";
      }
      ostream << "\n";
    }
    return ostream;
  }

  /**
   * \brief      Addition assignment operator with an expression
   *
   * \param      expr  The expression
   *
   * \tparam     E     The type of the expression
   *
   * \return     This matrix
   */
  template <typename E>
  inline matrix<T> &operator+=(expression<E> const &expr);

  /**
   * \brief      Multiplication assignment operator with an expression
   *
   * \param      expr  The expression
   *
   * \tparam     E     The type of the expression
   *
   * \return     This matrix
   */
  template <typename E>
  inline matrix<T> &operator*=(expression<E> const &expr);

  /**
   * \brief      Subtraction assignment operator with an expression
   *
   * \param      expr  The expression
   *
   * \tparam     E     The type of the expression
   *
   * \return     This matrix
   */
  template <typename E>
  inline matrix<T> &operator-=(expression<E> const &expr);

  /**
   * \brief      Addition assignment operator with a scalar
   *
   * \param      expr    The expression
   *
   * \tparam     Scalar  The type of the scalar
   *
   * \return     This matrix
   */
  template <typename Scalar, typename = enable_if_not_expression<Scalar>>
  inline matrix<T> &operator+=(Scalar const &expr);

  /**
   * \brief      Subtraction assignment operator with a scalar
   *
   * \param      expr    The expression
   *
   * \tparam     Scalar  The type of the scalar
   *
   * \return     This matrix
   */
  template <typename Scalar, typename = enable_if_not_expression<Scalar>>
  inline matrix<T> &operator-=(Scalar const &expr);

  /**
   * \brief      Multiplication assignment operator with a scalar
   *
   * \param      expr    The expression
   *
   * \tparam     Scalar  The type of the scalar
   *
   * \return     This matrix
   */
  template <typename Scalar, typename = enable_if_not_expression<Scalar>>
  inline matrix<T> &operator*=(Scalar const &expr);
};
} // namespace fastmatrix

namespace std {
using namespace fastmatrix;

/**
 * \brief      Overloading common_type trait for two matrices with different element types
 *
 * \tparam     T1    Type of elements of matrix 1
 * \tparam     T2    Type of elements of matrix 2
 */
template <typename T1, typename T2>
struct common_type<matrix<T1>, matrix<T2>> {
  using type = matrix<std::common_type_t<T1, T2>>;
};

/**
 * \brief      Overloading common_type trait for a matrix and a scalar
 *
 * \tparam     T1    Type of elements of matrix 1
 * \tparam     T2    Type of scalar
 */
template <typename T1, typename T2>
struct common_type<matrix<T1>, T2> {
  using type = matrix<std::common_type_t<T1, T2>>;
};
} // namespace std

namespace fastmatrix {

/**
 * \brief      Template specialization for storage type of scalar expressions
 *
 * \tparam     T     Type of scalar
 */
template <typename T>
struct storage_type<scalar_expression<T>> {
  using type = scalar_expression<T>;
};

/**
 * \brief      Class for coefficient-wise (element-wise) binary operations on matrix expressions
 *
 * \tparam     Op    Operation to be performed on the two expressions
 * \tparam     E1    Type of expression 1, must be matrix expression
 * \tparam     E2    Type of expression 2, can be matrix or scalar expression
 */
template <typename Op, typename E1, typename E2>
class cwise_matrix_binary_operation : public expression<cwise_matrix_binary_operation<Op, E1, E2>> {
private:
  /**
   * Expression 1
   */
  E1 const &expr1;

  /**
   * Expression 2, using storage_type trait as matrix expressions are stored as const references and
   * scalar expressions are stored by value
   */
  storage_type_t<E2> expr2;

public:
  /**
   * Return type of eval() method
   */
  using EvalReturnType = std::common_type_t<eval_return_type_t<E1>, eval_return_type_t<E2>>;

  /**
   * Type of elements of this expression
   */
  using ElementType = std::common_type_t<element_type_t<E1>, element_type_t<E2>>;

  /**
   * \brief      Constructor
   *
   * \param      expr1  Expression 1
   * \param      expr2  Expression 2
   */
  inline cwise_matrix_binary_operation(expression<E1> const &expr1, expression<E2> const &expr2)
      : expr1(expr1.get_const_derived()), expr2(expr2.get_const_derived()) {}

  /**
   * \brief      Function call operator to get an element of this matrix
   *
   * This uses lazy evaluation when possible to apply an operation only to the required elements of
   * expression 1 and 2
   *
   * \param[in]  i     Row number of the element to get
   * \param[in]  j     Column number of the element to get
   *
   * \return     The desired element
   */
  inline ElementType operator()(std::size_t i, std::size_t j) const {
    return Op::apply(expr1, expr2, i, j);
  }

  /**
   * \brief      Evaluates this expression fully and returns the resultant matrix
   *
   * This can also be triggered by a user of this library manually when they think that eagerly
   * evaluating an expression is beneficial, e.g.: a * (b + c) would benefit from eager evaluation
   * of b + c to prevent repeated calculations of the sum b + c. Thus, it could be re-written as
   * a * (b + c).eval()
   *
   * \return     The result of evaluating this expression
   */
  inline const EvalReturnType eval() const {
    EvalReturnType temp(num_rows(), num_cols());
    temp.assign((*this));
    return temp;
  }

  /**
   * \brief      Gets number of rows in this expression
   *
   * \return     Number of rows
   */
  inline std::size_t num_rows() const {
    return expr1.num_rows();
  }

  /**
   * \brief      Gets the number of columns in this expression
   *
   * \return     Number of columns
   */
  inline std::size_t num_cols() const {
    return expr1.num_cols();
  }
};

/**
 * \brief      Makes a cwise matrix binary operation
 *
 * Helps with template argument inference, so that the types of operands need not be specified while
 * specifying the operator
 *
 * \param      expr1  Expression 1
 * \param      expr2  Expression 2
 *
 * \tparam     Op     Type of operation
 * \tparam     E1     Type of expression 1
 * \tparam     E2     Type of expression 2
 *
 * \return     A cwise matrix binary operation
 */
template <template <typename E1, typename E2> typename Op, typename E1, typename E2>
cwise_matrix_binary_operation<Op<E1, E2>, E1, E2>
make_cwise_matrix_binary_operation(expression<E1> const &expr1, expression<E2> const &expr2) {
  return cwise_matrix_binary_operation<Op<E1, E2>, E1, E2>(expr1, expr2);
}

/**
 * \brief      Class to represent a matrix product
 *
 * Matrix products are eagerly evaluated and stored in a temporary, as statements like x = x * a can
 * cause incorrect results if lazily evaluated
 *
 * \tparam     E1    Type of expression 1
 * \tparam     E2    Type of expression 2
 */
template <typename E1, typename E2>
class matrix_product : public expression<matrix_product<E1, E2>> {
public:
  /**
   * Return type of eval() method
   */
  using EvalReturnType = std::common_type_t<eval_return_type_t<E1>, eval_return_type_t<E2>>;

  /**
   * Type of an element in this matrix product
   */
  using ElementType = std::common_type_t<element_type_t<E1>, element_type_t<E2>>;

private:
  /**
   * Expression 1
   */
  E1 const &expr1;

  /**
   * Expression 2
   */
  E2 const &expr2;

  /**
   * Temporary to store the product of the two expressions
   */
  EvalReturnType temp;

public:
  /**
   * \brief      Constructor
   *
   * Evaluates the matrix product of the two given expressions and stores it in a temporary
   *
   * \param      expr1  Expression 1
   * \param      expr2  Expression 2
   */
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

  /**
   * \brief      Function call operator to get an element of the matrix product
   *
   * \param[in]  i     Row number of the element to get
   * \param[in]  j     Column number of the element to get
   *
   * \return     The desired element
   */
  inline ElementType operator()(std::size_t i, std::size_t j) const {
    return temp(i, j);
  }

  /**
   * \brief      Evaluate and return the result of this expression
   *
   * Since the matrix product expression is eagerly evaluated as soon as it is constructed, this
   * method simply returns a const reference to the temporary in which the result has been stored
   *
   * \return     Reference to the product
   */
  inline EvalReturnType const &eval() const {
    return temp;
  }

  /**
   * \brief      Get number of rows of this matrix product
   *
   * \return     Number of rows
   */
  inline std::size_t num_rows() const {
    return expr1.num_rows();
  }

  /**
   * \brief      Get number of columns of this matrix product
   *
   * \return     Number of columns
   */
  inline std::size_t num_cols() const {
    return expr2.num_cols();
  }
};

/**
 * \brief      Operator representing element-wise matrix addition between two expressions
 *
 * \tparam     E1    Type of expression 1
 * \tparam     E2    Type of expression 2
 */
template <typename E1, typename E2>
struct cwise_matrix_add {
  /**
   * \brief      Applies the addition operator on the given elements of the given expressions
   *
   * \param      expr1  Expression 1
   * \param      expr2  Expression 2
   * \param[in]  i      Row number of the element
   * \param[in]  j      Column number of the element
   *
   * \return     Sum of the two elements
   */
  inline static typename std::common_type_t<element_type_t<E1>, element_type_t<E2>>
  apply(expression<E1> const &expr1, expression<E2> const &expr2, std::size_t i, std::size_t j) {
    return expr1.get_const_derived()(i, j) + expr2.get_const_derived()(i, j);
  }
};

/**
 * \brief      Operator representing element-wise matrix multiplication between two expressions
 *
 * \tparam     E1    Type of expression 1
 * \tparam     E2    Type of expression 2
 */
template <typename E1, typename E2>
struct cwise_matrix_multiply {
  /**
   * \brief      Applies the coefficient-wise multiplication operator on the given elements of the
   * given expressions
   *
   * \param      expr1  Expression 1
   * \param      expr2  Expression 1
   * \param[in]  i      Row number of the element
   * \param[in]  j      Column number of the element
   *
   * \return     Product of the two elements
   */
  inline static typename std::common_type_t<element_type_t<E1>, element_type_t<E2>>
  apply(expression<E1> const &expr1, expression<E2> const &expr2, std::size_t i, std::size_t j) {
    return expr1.get_const_derived()(i, j) * expr2.get_const_derived()(i, j);
  }
};

/**
 * \brief      Operator representing element-wise matrix subtraction between two expressions
 *
 * \tparam     E1    Type of expression 1
 * \tparam     E2    Type of expression 2
 */
template <typename E1, typename E2>
struct cwise_matrix_subtract {
  /**
   * \brief      Applies the subtraction operator on the given elements of the given expressions
   *
   * \param      expr1  Expression 1
   * \param      expr2  Expression 2
   * \param[in]  i      Row number of the element
   * \param[in]  j      Column number of the element
   *
   * \return     Difference of the two elements
   */
  inline static typename std::common_type_t<element_type_t<E1>, element_type_t<E2>>
  apply(expression<E1> const &expr1, expression<E2> const &expr2, std::size_t i, std::size_t j) {
    return expr1.get_const_derived()(i, j) - expr2.get_const_derived()(i, j);
  }
};

// The below functions and methods are simply arithmetic operator overloads to easily construct
// expression classes. e.g. a + b constructs a cwise_matrix_binary_operation with cwise_matrix_add
// as the operation

template <typename E1, typename E2>
inline auto operator+(expression<E1> const &expr1, expression<E2> const &expr2) {
  assert(expr1.num_rows() == expr2.num_rows());
  assert(expr1.num_cols() == expr2.num_cols());
  return make_cwise_matrix_binary_operation<cwise_matrix_add>(expr1, expr2);
}

template <typename E, typename T, typename = enable_if_not_expression<T>>
inline auto operator+(expression<E> const &expr, T const &scalar) {
  return make_cwise_matrix_binary_operation<cwise_matrix_add>(expr, scalar_expression(scalar));
}

template <typename E, typename T, typename = enable_if_not_expression<T>>
inline auto operator+(T const &scalar, expression<E> const &expr) {
  return expr + scalar;
}

template <typename E1, typename E2>
inline auto operator*(expression<E1> const &expr1, expression<E2> const &expr2) {
  assert(expr1.num_cols() == expr2.num_rows());
  return matrix_product(expr1, expr2);
}

template <typename E, typename T, typename = enable_if_not_expression<T>>
inline auto operator*(expression<E> const &expr, T const &scalar) {
  return make_cwise_matrix_binary_operation<cwise_matrix_multiply>(expr, scalar_expression(scalar));
}

template <typename E, typename T, typename = enable_if_not_expression<T>>
inline auto operator*(T const &scalar, expression<E> const &expr) {
  return expr * scalar;
}

template <typename E1, typename E2>
inline auto operator-(expression<E1> const &expr1, expression<E2> const &expr2) {
  assert(expr1.num_rows() == expr2.num_rows());
  assert(expr1.num_cols() == expr2.num_cols());
  return make_cwise_matrix_binary_operation<cwise_matrix_subtract>(expr1, expr2);
}

template <typename E, typename T, typename = enable_if_not_expression<T>>
inline auto operator-(expression<E> const &expr, T const &scalar) {
  return make_cwise_matrix_binary_operation<cwise_matrix_subtract>(expr, scalar_expression(scalar));
}

template <typename E, typename T, typename = enable_if_not_expression<T>>
inline auto operator-(T const &scalar, expression<E> const &expr) {
  return expr - scalar;
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

#endif // FASTMATRIX_FASTMATRIX_HPP
