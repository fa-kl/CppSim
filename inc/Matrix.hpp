/*****************************************************************************************
 * @file: Matrix.hpp
 *
 * @brief: This file provides a mathematical matrix class.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include <ostream>
#include <vector>

#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

/**
 * @brief A class for real-valued matrices.
 *
 * @details This class provides a mathematical matrix implementation with operations
 * for linear algebra, including matrix multiplication, decompositions, and element-wise functions.
 */
class Matrix
{
protected:
  /**
   * @brief The underlying data storage (row-major order).
   */
  std::vector<real_t> m_data;

  /**
   * @brief The dimensions of the matrix.
   */
  dimension_t m_size;

public:
  /**
   * @brief Creates an empty matrix.
   */
  Matrix();

  /**
   * @brief Creates a matrix with specified dimensions initialized to zero.
   *
   * @param rows Number of rows.
   * @param cols Number of columns.
   */
  Matrix(size_t rows, size_t cols);

  /**
   * @brief Creates a matrix from a nested initializer list.
   *
   * @param values Nested initializer list of values.
   */
  Matrix(std::initializer_list<std::initializer_list<real_t>> values);

  /**
   * @brief Creates a matrix from a vector by reshaping it.
   *
   * @param vec Vector to convert to matrix.
   * @param rows Number of rows.
   * @param cols Number of columns.
   */
  Matrix(const Vector& vec, size_t rows, size_t cols);

  /**
   * @brief Copy constructor.
   *
   * @param other Matrix to copy.
   */
  Matrix(const Matrix& other);

  /**
   * @brief Move constructor.
   *
   * @param other Matrix to move.
   */
  Matrix(Matrix&& other) noexcept;

  /**
   * @brief Destructor.
   */
  ~Matrix();

  /**
   * @brief Copy assignment operator.
   *
   * @param other Matrix to copy.
   * @return Reference to this matrix.
   */
  Matrix& operator=(const Matrix& other);

  /**
   * @brief Move assignment operator.
   *
   * @param other Matrix to move.
   * @return Reference to this matrix.
   */
  Matrix& operator=(Matrix&& other) noexcept;

  /**
   * @brief Get the number of rows.
   *
   * @return Number of rows.
   */
  size_t rows() const;

  /**
   * @brief Get the number of columns.
   *
   * @return Number of columns.
   */
  size_t cols() const;

  /**
   * @brief Get the size dimensions of the matrix.
   *
   * @return Size as dimension_t structure.
   */
  dimension_t size() const;

  /**
   * @brief Get the total number of elements.
   *
   * @return Total number of elements.
   */
  size_t length() const;

  /**
   * @brief Check if the matrix is empty.
   *
   * @return True if empty, false otherwise.
   */
  bool_t isEmpty() const;

  /**
   * @brief Access element at linear index (0-based, row-major order).
   *
   * @param i Linear index.
   * @return Reference to element at index i.
   */
  real_t& operator[](size_t i);

  /**
   * @brief Access element at linear index (0-based, row-major order).
   *
   * @param i Linear index.
   * @return Const reference to element at index i.
   */
  const real_t& operator[](size_t i) const;

  /**
   * @brief Access element at row i, column j (1-based indexing, supports negative indices).
   *
   * @param i Row index (1-based, negative indices count from end).
   * @param j Column index (1-based, negative indices count from end).
   * @return Reference to element at (i, j).
   */
  real_t& operator()(index_t i, index_t j);

  /**
   * @brief Access element at row i, column j (1-based indexing, supports negative indices).
   *
   * @param i Row index (1-based, negative indices count from end).
   * @param j Column index (1-based, negative indices count from end).
   * @return Const reference to element at (i, j).
   */
  const real_t& operator()(index_t i, index_t j) const;

  /**
   * @brief Unary plus operator.
   *
   * @return A copy of the matrix.
   */
  Matrix operator+() const;

  /**
   * @brief Unary minus operator (negate matrix).
   *
   * @return A negated copy of the matrix.
   */
  Matrix operator-() const;

  /**
   * @brief Add two matrices element-wise.
   *
   * @param rhs Right-hand side matrix.
   * @return Result of element-wise addition.
   */
  Matrix operator+(const Matrix& rhs) const;

  /**
   * @brief Add a scalar to all elements of the matrix.
   *
   * @param scalar Scalar value.
   * @return Result matrix.
   */
  Matrix operator+(real_t scalar) const;

  /**
   * @brief Subtract two matrices element-wise.
   *
   * @param rhs Right-hand side matrix.
   * @return Result of element-wise subtraction.
   */
  Matrix operator-(const Matrix& rhs) const;

  /**
   * @brief Subtract a scalar from all elements of the matrix.
   *
   * @param scalar Scalar value.
   * @return Result matrix.
   */
  Matrix operator-(real_t scalar) const;

  /**
   * @brief Perform matrix multiplication.
   *
   * @param rhs Right-hand side matrix.
   * @return Result of matrix multiplication.
   */
  Matrix operator*(const Matrix& rhs) const;

  /**
   * @brief Multiply matrix by a vector.
   *
   * @param vec Vector.
   * @return Result vector.
   */
  Vector operator*(const Vector& vec) const;

  /**
   * @brief Multiply all elements by a scalar.
   *
   * @param scalar Scalar value.
   * @return Result matrix.
   */
  Matrix operator*(real_t scalar) const;

  /**
   * @brief Divide all elements by a scalar.
   *
   * @param scalar Scalar value (divisor).
   * @return Result matrix.
   */
  Matrix operator/(real_t scalar) const;

  /**
   * @brief Element-wise equality comparison.
   *
   * @param rhs Right-hand side matrix.
   * @return True if all elements are equal, false otherwise.
   */
  bool_t operator==(const Matrix& rhs) const;

  /**
   * @brief Element-wise inequality comparison.
   *
   * @param rhs Right-hand side matrix.
   * @return True if any elements differ, false otherwise.
   */
  bool_t operator!=(const Matrix& rhs) const;

  /**
   * @brief Multiply a scalar by a matrix.
   *
   * @param scalar Scalar value.
   * @param mat Matrix.
   * @return Result matrix.
   */
  friend Matrix operator*(real_t scalar, const Matrix& mat);

  /**
   * @brief Multiply a vector by a matrix.
   *
   * @param vec Vector.
   * @param mat Matrix.
   * @return Result matrix (outer product-like operation).
   */
  friend Matrix operator*(const Vector& vec, const Matrix& mat);

  /**
   * @brief Add a scalar to a matrix.
   *
   * @param scalar Scalar value.
   * @param mat Matrix.
   * @return Result matrix.
   */
  friend Matrix operator+(real_t scalar, const Matrix& mat);

  /**
   * @brief Subtract a matrix from a scalar.
   *
   * @param scalar Scalar value.
   * @param mat Matrix.
   * @return Result matrix.
   */
  friend Matrix operator-(real_t scalar, const Matrix& mat);

  /**
   * @brief Output stream operator for matrices.
   *
   * @param os Output stream.
   * @param mat Matrix to output.
   * @return Output stream.
   */
  friend std::ostream& operator<<(std::ostream& os, const Matrix& mat);
};

/*****************************************************************************************
 * Matrix Functions
 ****************************************************************************************/

/**
 * @brief Create a matrix of zeros.
 *
 * @param rows Number of rows.
 * @param cols Number of columns.
 * @return Matrix filled with zeros.
 */
Matrix zeros(size_t rows, size_t cols);

/**
 * @brief Create a matrix of ones.
 *
 * @param rows Number of rows.
 * @param cols Number of columns.
 * @return Matrix filled with ones.
 */
Matrix ones(size_t rows, size_t cols);

/**
 * @brief Create an identity matrix.
 *
 * @param n Size of the square matrix.
 * @return Identity matrix of size n x n.
 */
Matrix eye(size_t n);

/**
 * @brief Create a diagonal matrix from an initializer list.
 *
 * @param values Diagonal values.
 * @return Diagonal matrix.
 */
Matrix diag(std::initializer_list<real_t> values);

/**
 * @brief Create a diagonal matrix from a vector.
 *
 * @param vec Vector containing diagonal values.
 * @return Diagonal matrix.
 */
Matrix diag(const Vector& vec);

/**
 * @brief Extract the diagonal of a matrix as a vector.
 *
 * @param mat Input matrix.
 * @return Vector containing diagonal elements.
 */
Vector diag(const Matrix& mat);

/**
 * @brief Get the number of rows in a matrix.
 *
 * @param mat Input matrix.
 * @return Number of rows.
 */
size_t rows(const Matrix& mat);

/**
 * @brief Get the number of columns in a matrix.
 *
 * @param mat Input matrix.
 * @return Number of columns.
 */
size_t cols(const Matrix& mat);

/**
 * @brief Get the size dimensions of a matrix.
 *
 * @param mat Input matrix.
 * @return Size as dimension_t structure.
 */
dimension_t size(const Matrix& mat);

/**
 * @brief Get the total number of elements in a matrix.
 *
 * @param mat Input matrix.
 * @return Total number of elements.
 */
size_t length(const Matrix& mat);

/**
 * @brief Check if a matrix is empty.
 *
 * @param mat Input matrix.
 * @return True if empty, false otherwise.
 */
bool_t isEmpty(const Matrix& mat);

/**
 * @brief Transpose a matrix (swap rows and columns).
 *
 * @param mat Input matrix.
 * @return Transposed matrix.
 */
Matrix transpose(const Matrix& mat);

/**
 * @brief Transpose a vector (convert to row matrix).
 *
 * @param vec Input vector.
 * @return Transposed matrix (1 x n).
 */
Matrix transpose(const Vector& vec);

/**
 * @brief Compute the trace of a matrix (sum of diagonal elements).
 *
 * @param mat Input matrix.
 * @return Trace value.
 */
real_t trace(const Matrix& mat);

/**
 * @brief Compute absolute value of each element in a matrix.
 *
 * @param mat Input matrix.
 * @return Matrix with absolute values.
 */
Matrix abs(const Matrix& mat);

/**
 * @brief Compute exponential of each element in a matrix.
 *
 * @param mat Input matrix.
 * @return Matrix with exponential values.
 */
Matrix exp(const Matrix& mat);

/**
 * @brief Compute natural logarithm of each element in a matrix.
 *
 * @param mat Input matrix.
 * @return Matrix with logarithmic values.
 */
Matrix log(const Matrix& mat);

/**
 * @brief Compute sine of each element in a matrix.
 *
 * @param mat Input matrix.
 * @return Matrix with sine values.
 */
Matrix sin(const Matrix& mat);

/**
 * @brief Compute cosine of each element in a matrix.
 *
 * @param mat Input matrix.
 * @return Matrix with cosine values.
 */
Matrix cos(const Matrix& mat);

/**
 * @brief Compute tangent of each element in a matrix.
 *
 * @param mat Input matrix.
 * @return Matrix with tangent values.
 */
Matrix tan(const Matrix& mat);

/**
 * @brief Compute arcsine of each element in a matrix.
 *
 * @param mat Input matrix.
 * @return Matrix with arcsine values.
 */
Matrix asin(const Matrix& mat);

/**
 * @brief Compute arccosine of each element in a matrix.
 *
 * @param mat Input matrix.
 * @return Matrix with arccosine values.
 */
Matrix acos(const Matrix& mat);

/**
 * @brief Compute arctangent of each element in a matrix.
 *
 * @param mat Input matrix.
 * @return Matrix with arctangent values.
 */
Matrix atan(const Matrix& mat);

/**
 * @brief Raise each element of a matrix to a power (element-wise).
 *
 * @param mat Input matrix.
 * @param exponent Power exponent.
 * @return Matrix with powered values.
 */
Matrix elpow(const Matrix& mat, real_t exponent);

/**
 * @brief Compute square root of each element in a matrix.
 *
 * @param mat Input matrix.
 * @return Matrix with square root values.
 */
Matrix sqrt(const Matrix& mat);

/**
 * @brief Compute sum of all elements in a matrix.
 *
 * @param mat Input matrix.
 * @return Sum of all elements.
 */
real_t sum(const Matrix& mat);

/**
 * @brief Element-wise multiplication of two matrices.
 *
 * @param lhs Left-hand side matrix.
 * @param rhs Right-hand side matrix.
 * @return Matrix with element-wise products.
 */
Matrix elmul(const Matrix& lhs, const Matrix& rhs);

/**
 * @brief Element-wise division of two matrices.
 *
 * @param lhs Left-hand side matrix (dividend).
 * @param rhs Right-hand side matrix (divisor).
 * @return Matrix with element-wise quotients.
 */
Matrix eldiv(const Matrix& lhs, const Matrix& rhs);

/**
 * @brief Find maximum element in a matrix.
 *
 * @param mat Input matrix.
 * @return Maximum value.
 */
real_t max(const Matrix& mat);

/**
 * @brief Find minimum element in a matrix.
 *
 * @param mat Input matrix.
 * @return Minimum value.
 */
real_t min(const Matrix& mat);

/**
 * @brief Compute mean (average) of all matrix elements.
 *
 * @param mat Input matrix.
 * @return Mean value.
 */
real_t mean(const Matrix& mat);

/**
 * @brief Compute variance of all matrix elements.
 *
 * @param mat Input matrix.
 * @return Variance.
 */
real_t var(const Matrix& mat);

/**
 * @brief Compute standard deviation of all matrix elements.
 *
 * @param mat Input matrix.
 * @return Standard deviation.
 */
real_t std(const Matrix& mat);

/**
 * @brief Horizontally concatenate two matrices.
 *
 * @param lhs Left matrix.
 * @param rhs Right matrix.
 * @return Concatenated matrix.
 */
Matrix hcat(const Matrix& lhs, const Matrix& rhs);

/**
 * @brief Horizontally concatenate two vectors into a matrix.
 *
 * @param lhs Left vector (first column).
 * @param rhs Right vector (second column).
 * @return Matrix with two columns.
 */
Matrix hcat(const Vector& lhs, const Vector& rhs);

/**
 * @brief Horizontally concatenate a matrix and a vector.
 *
 * @param lhs Matrix.
 * @param rhs Vector.
 * @return Concatenated matrix.
 */
Matrix hcat(const Matrix& lhs, const Vector& rhs);

/**
 * @brief Horizontally concatenate a vector and a matrix.
 *
 * @param lhs Vector.
 * @param rhs Matrix.
 * @return Concatenated matrix.
 */
Matrix hcat(const Vector& lhs, const Matrix& rhs);

/**
 * @brief Vertically concatenate two matrices.
 *
 * @param lhs Top matrix.
 * @param rhs Bottom matrix.
 * @return Concatenated matrix.
 */
Matrix vcat(const Matrix& lhs, const Matrix& rhs);

/**
 * @brief Compute the Frobenius norm of a matrix.
 *
 * @param mat Input matrix.
 * @param p Norm parameter (default: 2).
 * @return Frobenius norm.
 */
real_t norm(const Matrix& mat, real_t p = 2.0);

/**
 * @brief Compute the outer product of two vectors.
 *
 * @param lhs Left-hand side vector.
 * @param rhs Right-hand side vector.
 * @return Matrix representing the outer product.
 */
Matrix outer(const Vector& lhs, const Vector& rhs);

/**
 * @brief Compute the inverse of a square matrix using Gaussian elimination.
 *
 * @param mat Input square matrix.
 * @return Inverse matrix.
 */
Matrix inv(const Matrix& mat);

}  // namespace sim