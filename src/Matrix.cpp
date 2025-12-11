/*****************************************************************************************
 * @file: Matrix.cpp
 *
 * @brief: This file implements the matrix's functions.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Matrix.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Vector.hpp"

namespace sim
{

/*****************************************************************************************
 * Constructors and Destructors
 ****************************************************************************************/

Matrix::Matrix() : m_data(), m_size({0, 0}) {}

Matrix::Matrix(size_t rows, size_t cols) : m_data(rows * cols, 0.0), m_size({rows, cols}) {}

Matrix::Matrix(std::initializer_list<std::initializer_list<real_t>> values)
    : m_data(), m_size({values.size(), values.size() > 0 ? values.begin()->size() : 0})
{
  m_data.reserve(m_size.rows * m_size.cols);
  for (const auto& row : values) {
    if (row.size() != m_size.cols) {
      throw std::invalid_argument("All rows must have the same number of columns");
    }
    m_data.insert(m_data.end(), row.begin(), row.end());
  }
}

Matrix::Matrix(const Vector& vec, size_t rows, size_t cols) : m_data(rows * cols), m_size({rows, cols})
{
  if (vec.length() != rows * cols) {
    throw std::invalid_argument("Vector length must match matrix dimensions");
  }
  for (size_t i = 0; i < vec.length(); ++i) {
    m_data[i] = vec[i];
  }
}
Matrix::Matrix(const Matrix& other) : m_data(other.m_data), m_size(other.m_size) {}

Matrix::Matrix(Matrix&& other) noexcept : m_data(std::move(other.m_data)), m_size(other.m_size)
{
  other.m_size = {0, 0};
}

Matrix::~Matrix() = default;

/*****************************************************************************************
 * Assignment Operators
 ****************************************************************************************/

Matrix& Matrix::operator=(const Matrix& other)
{
  if (this != &other) {
    m_data = other.m_data;
    m_size = other.m_size;
  }
  return *this;
}

Matrix& Matrix::operator=(Matrix&& other) noexcept
{
  if (this != &other) {
    m_data = std::move(other.m_data);
    m_size = other.m_size;
    other.m_size = {0, 0};
  }
  return *this;
}

/*****************************************************************************************
 * Size and Capacity
 ****************************************************************************************/

size_t Matrix::rows() const
{
  return m_size.rows;
}

size_t Matrix::cols() const
{
  return m_size.cols;
}

dimension_t Matrix::size() const
{
  return m_size;
}

size_t Matrix::length() const
{
  return m_size.rows * m_size.cols;
}

bool_t Matrix::isEmpty() const
{
  return m_data.empty();
}

/*****************************************************************************************
 * Element Access
 ****************************************************************************************/

real_t& Matrix::operator[](size_t i)
{
  return m_data[i];
}

const real_t& Matrix::operator[](size_t i) const
{
  return m_data[i];
}

real_t& Matrix::operator()(index_t i, index_t j)
{
  if (i == 0 || j == 0) {
    throw std::out_of_range("Indices cannot be zero (1-based indexing)");
  }
  size_t row = (i > 0) ? static_cast<size_t>(i - 1) : (m_size.rows + static_cast<size_t>(i));
  size_t col = (j > 0) ? static_cast<size_t>(j - 1) : (m_size.cols + static_cast<size_t>(j));
  if (row >= m_size.rows || col >= m_size.cols) {
    throw std::out_of_range("Index out of range");
  }
  return m_data[row * m_size.cols + col];
}

const real_t& Matrix::operator()(index_t i, index_t j) const
{
  if (i == 0 || j == 0) {
    throw std::out_of_range("Indices cannot be zero (1-based indexing)");
  }
  size_t row = (i > 0) ? static_cast<size_t>(i - 1) : (m_size.rows + static_cast<size_t>(i));
  size_t col = (j > 0) ? static_cast<size_t>(j - 1) : (m_size.cols + static_cast<size_t>(j));
  if (row >= m_size.rows || col >= m_size.cols) {
    throw std::out_of_range("Index out of range");
  }
  return m_data[row * m_size.cols + col];
}

/*****************************************************************************************
 * Unary Operators
 ****************************************************************************************/

Matrix Matrix::operator+() const
{
  return *this;
}

Matrix Matrix::operator-() const
{
  Matrix result(m_size.rows, m_size.cols);
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = -m_data[i];
  }
  return result;
}

/*****************************************************************************************
 * Arithmetic Operators
 ****************************************************************************************/

Matrix Matrix::operator+(const Matrix& rhs) const
{
  if (m_size.rows != rhs.m_size.rows || m_size.cols != rhs.m_size.cols) {
    throw std::invalid_argument("Incompatible matrix sizes");
  }
  Matrix result(m_size.rows, m_size.cols);
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = m_data[i] + rhs.m_data[i];
  }
  return result;
}

Matrix Matrix::operator+(real_t scalar) const
{
  Matrix result(m_size.rows, m_size.cols);
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = m_data[i] + scalar;
  }
  return result;
}

Matrix Matrix::operator-(const Matrix& rhs) const
{
  if (m_size.rows != rhs.m_size.rows || m_size.cols != rhs.m_size.cols) {
    throw std::invalid_argument("Incompatible matrix sizes");
  }
  Matrix result(m_size.rows, m_size.cols);
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = m_data[i] - rhs.m_data[i];
  }
  return result;
}

Matrix Matrix::operator-(real_t scalar) const
{
  Matrix result(m_size.rows, m_size.cols);
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = m_data[i] - scalar;
  }
  return result;
}

Matrix Matrix::operator*(const Matrix& rhs) const
{
  if (m_size.cols != rhs.m_size.rows) {
    throw std::invalid_argument("Incompatible matrix dimensions for multiplication");
  }
  Matrix result(m_size.rows, rhs.m_size.cols);
  for (size_t i = 0; i < m_size.rows; ++i) {
    for (size_t j = 0; j < rhs.m_size.cols; ++j) {
      real_t sum = 0.0;
      for (size_t k = 0; k < m_size.cols; ++k) {
        sum += m_data[i * m_size.cols + k] * rhs.m_data[k * rhs.m_size.cols + j];
      }
      result.m_data[i * rhs.m_size.cols + j] = sum;
    }
  }
  return result;
}

Vector Matrix::operator*(const Vector& vec) const
{
  if (m_size.cols != vec.length()) {
    throw std::invalid_argument("Incompatible dimensions for matrix-vector multiplication");
  }
  Vector result(m_size.rows);
  for (size_t i = 0; i < m_size.rows; ++i) {
    real_t sum = 0.0;
    for (size_t j = 0; j < m_size.cols; ++j) {
      sum += m_data[i * m_size.cols + j] * vec[j];
    }
    result[i] = sum;
  }
  return result;
}

Matrix Matrix::operator*(real_t scalar) const
{
  Matrix result(m_size.rows, m_size.cols);
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = m_data[i] * scalar;
  }
  return result;
}

Matrix Matrix::operator/(real_t scalar) const
{
  if (std::abs(scalar) < 1e-10) {
    throw std::invalid_argument("Division by zero");
  }
  Matrix result(m_size.rows, m_size.cols);
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = m_data[i] / scalar;
  }
  return result;
}

Matrix operator*(real_t scalar, const Matrix& mat)
{
  return mat * scalar;
}

Matrix operator*(const Vector& vec, const Matrix& mat)
{
  if (vec.length() != mat.m_size.rows) {
    throw std::invalid_argument("Incompatible dimensions for vector-matrix multiplication");
  }
  Matrix result(1, mat.m_size.cols);
  for (size_t j = 0; j < mat.m_size.cols; ++j) {
    real_t sum = 0.0;
    for (size_t i = 0; i < vec.length(); ++i) {
      sum += vec[i] * mat.m_data[i * mat.m_size.cols + j];
    }
    result.m_data[j] = sum;
  }
  return result;
}

Matrix operator+(real_t scalar, const Matrix& mat)
{
  return mat + scalar;
}

Matrix operator-(real_t scalar, const Matrix& mat)
{
  Matrix result(mat.m_size.rows, mat.m_size.cols);
  for (size_t i = 0; i < mat.m_data.size(); ++i) {
    result.m_data[i] = scalar - mat.m_data[i];
  }
  return result;
}

/*****************************************************************************************
 * Comparison Operators
 ****************************************************************************************/

bool_t Matrix::operator==(const Matrix& rhs) const
{
  if (m_size.rows != rhs.m_size.rows || m_size.cols != rhs.m_size.cols) {
    return false;
  }
  for (size_t i = 0; i < m_data.size(); ++i) {
    if (m_data[i] != rhs.m_data[i]) {
      return false;
    }
  }
  return true;
}

bool_t Matrix::operator!=(const Matrix& rhs) const
{
  return !(*this == rhs);
}

/*****************************************************************************************
 * Print Function(s)
 ****************************************************************************************/

std::ostream& operator<<(std::ostream& os, const Matrix& mat)
{
  // Prepare formatted strings with fixed 3 decimal places per element
  size_t rows = mat.rows();
  size_t cols = mat.cols();
  std::vector<std::vector<std::string>> parts(rows, std::vector<std::string>(cols));
  std::vector<size_t> col_width(cols, 0);

  std::ostringstream tmp;
  tmp.setf(std::ios::fixed, std::ios::floatfield);
  tmp.precision(3);

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      tmp.str("");
      tmp.clear();
      tmp << mat(static_cast<index_t>(i + 1), static_cast<index_t>(j + 1));
      parts[i][j] = tmp.str();
      col_width[j] = std::max(col_width[j], parts[i][j].size());
    }
  }

  os << "[";
  for (size_t i = 0; i < rows; ++i) {
    if (i > 0)
      os << " ";
    for (size_t j = 0; j < cols; ++j) {
      os << std::setw(static_cast<int>(col_width[j])) << parts[i][j];
      if (j + 1 < cols)
        os << ", ";
    }
    if (i + 1 < rows)
      os << "\n";
  }
  os << "]";
  return os;
}

/*******************************************************************************************
 * Matrix Functions
 ******************************************************************************************/

Matrix zeros(size_t rows, size_t cols)
{
  return Matrix(rows, cols);
}

Matrix ones(size_t rows, size_t cols)
{
  Matrix result(rows, cols);
  for (size_t i = 0; i < result.length(); ++i) {
    result[i] = 1.0;
  }
  return result;
}

Matrix eye(size_t n)
{
  Matrix result(n, n);
  for (size_t i = 0; i < n; ++i) {
    result[i * n + i] = 1.0;
  }
  return result;
}

Matrix diag(std::initializer_list<real_t> values)
{
  size_t n = values.size();
  Matrix result(n, n);
  size_t i = 0;
  for (real_t val : values) {
    result[i * n + i] = val;
    ++i;
  }
  return result;
}

Matrix diag(const Vector& vec)
{
  size_t n = vec.length();
  Matrix result(n, n);
  for (size_t i = 0; i < n; ++i) {
    result[i * n + i] = vec[i];
  }
  return result;
}

Vector diag(const Matrix& mat)
{
  size_t n = std::min(mat.rows(), mat.cols());
  Vector result(n);
  for (size_t i = 0; i < n; ++i) {
    result[i] = mat[i * mat.cols() + i];
  }
  return result;
}

size_t rows(const Matrix& mat)
{
  return mat.rows();
}

size_t cols(const Matrix& mat)
{
  return mat.cols();
}

dimension_t size(const Matrix& mat)
{
  return mat.size();
}

size_t length(const Matrix& mat)
{
  return mat.length();
}

bool_t isEmpty(const Matrix& mat)
{
  return mat.isEmpty();
}

Matrix transpose(const Matrix& mat)
{
  Matrix result(mat.cols(), mat.rows());
  for (size_t i = 0; i < mat.rows(); ++i) {
    for (size_t j = 0; j < mat.cols(); ++j) {
      result[j * mat.rows() + i] = mat[i * mat.cols() + j];
    }
  }
  return result;
}

Matrix transpose(const Vector& vec)
{
  Matrix result(1, vec.length());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = vec[i];
  }
  return result;
}

real_t trace(const Matrix& mat)
{
  size_t n = std::min(mat.rows(), mat.cols());
  real_t sum = 0.0;
  for (size_t i = 0; i < n; ++i) {
    sum += mat[i * mat.cols() + i];
  }
  return sum;
}

Matrix abs(const Matrix& mat)
{
  Matrix result(mat.rows(), mat.cols());
  for (size_t i = 0; i < mat.length(); ++i) {
    result[i] = std::abs(mat[i]);
  }
  return result;
}

Matrix exp(const Matrix& mat)
{
  Matrix result(mat.rows(), mat.cols());
  for (size_t i = 0; i < mat.length(); ++i) {
    result[i] = std::exp(mat[i]);
  }
  return result;
}

Matrix log(const Matrix& mat)
{
  Matrix result(mat.rows(), mat.cols());
  for (size_t i = 0; i < mat.length(); ++i) {
    result[i] = std::log(mat[i]);
  }
  return result;
}

Matrix sin(const Matrix& mat)
{
  Matrix result(mat.rows(), mat.cols());
  for (size_t i = 0; i < mat.length(); ++i) {
    result[i] = std::sin(mat[i]);
  }
  return result;
}

Matrix cos(const Matrix& mat)
{
  Matrix result(mat.rows(), mat.cols());
  for (size_t i = 0; i < mat.length(); ++i) {
    result[i] = std::cos(mat[i]);
  }
  return result;
}

Matrix tan(const Matrix& mat)
{
  Matrix result(mat.rows(), mat.cols());
  for (size_t i = 0; i < mat.length(); ++i) {
    result[i] = std::tan(mat[i]);
  }
  return result;
}

Matrix asin(const Matrix& mat)
{
  Matrix result(mat.rows(), mat.cols());
  for (size_t i = 0; i < mat.length(); ++i) {
    result[i] = std::asin(mat[i]);
  }
  return result;
}

Matrix acos(const Matrix& mat)
{
  Matrix result(mat.rows(), mat.cols());
  for (size_t i = 0; i < mat.length(); ++i) {
    result[i] = std::acos(mat[i]);
  }
  return result;
}

Matrix atan(const Matrix& mat)
{
  Matrix result(mat.rows(), mat.cols());
  for (size_t i = 0; i < mat.length(); ++i) {
    result[i] = std::atan(mat[i]);
  }
  return result;
}

Matrix elpow(const Matrix& mat, real_t exponent)
{
  Matrix result(mat.rows(), mat.cols());
  for (size_t i = 0; i < mat.length(); ++i) {
    result[i] = std::pow(mat[i], exponent);
  }
  return result;
}

Matrix sqrt(const Matrix& mat)
{
  Matrix result(mat.rows(), mat.cols());
  for (size_t i = 0; i < mat.length(); ++i) {
    result[i] = std::sqrt(mat[i]);
  }
  return result;
}

real_t sum(const Matrix& mat)
{
  real_t result = 0.0;
  for (size_t i = 0; i < mat.length(); ++i) {
    result += mat[i];
  }
  return result;
}

Matrix elmul(const Matrix& lhs, const Matrix& rhs)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw std::invalid_argument("Incompatible matrix sizes");
  }
  Matrix result(lhs.rows(), lhs.cols());
  for (size_t i = 0; i < lhs.length(); ++i) {
    result[i] = lhs[i] * rhs[i];
  }
  return result;
}

Matrix eldiv(const Matrix& lhs, const Matrix& rhs)
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
    throw std::invalid_argument("Incompatible matrix sizes");
  }
  Matrix result(lhs.rows(), lhs.cols());
  for (size_t i = 0; i < lhs.length(); ++i) {
    if (std::abs(rhs[i]) < 1e-10) {
      throw std::invalid_argument("Division by zero");
    }
    result[i] = lhs[i] / rhs[i];
  }
  return result;
}

real_t max(const Matrix& mat)
{
  if (mat.isEmpty()) {
    throw std::invalid_argument("Cannot find max of empty matrix");
  }
  real_t max_val = mat[0];
  for (size_t i = 1; i < mat.length(); ++i) {
    max_val = std::max(max_val, mat[i]);
  }
  return max_val;
}

real_t min(const Matrix& mat)
{
  if (mat.isEmpty()) {
    throw std::invalid_argument("Cannot find min of empty matrix");
  }
  real_t min_val = mat[0];
  for (size_t i = 1; i < mat.length(); ++i) {
    min_val = std::min(min_val, mat[i]);
  }
  return min_val;
}

real_t mean(const Matrix& mat)
{
  if (mat.isEmpty()) {
    return 0.0;
  }
  return sum(mat) / static_cast<real_t>(mat.length());
}

real_t var(const Matrix& mat)
{
  if (mat.isEmpty()) {
    return 0.0;
  }
  real_t m = mean(mat);
  real_t sum_sq = 0.0;
  for (size_t i = 0; i < mat.length(); ++i) {
    real_t diff = mat[i] - m;
    sum_sq += diff * diff;
  }
  return sum_sq / static_cast<real_t>(mat.length());
}

real_t std(const Matrix& mat)
{
  return std::sqrt(var(mat));
}

Matrix hcat(const Matrix& lhs, const Matrix& rhs)
{
  if (lhs.rows() != rhs.rows()) {
    throw std::invalid_argument("Matrices must have the same number of rows");
  }
  Matrix result(lhs.rows(), lhs.cols() + rhs.cols());
  for (size_t i = 0; i < lhs.rows(); ++i) {
    for (size_t j = 0; j < lhs.cols(); ++j) {
      result[i * result.cols() + j] = lhs[i * lhs.cols() + j];
    }
    for (size_t j = 0; j < rhs.cols(); ++j) {
      result[i * result.cols() + lhs.cols() + j] = rhs[i * rhs.cols() + j];
    }
  }
  return result;
}

Matrix hcat(const Vector& lhs, const Vector& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw std::invalid_argument("Vectors must have the same length");
  }
  Matrix result(lhs.length(), 2);
  for (size_t i = 0; i < lhs.length(); ++i) {
    result[i * 2] = lhs[i];
    result[i * 2 + 1] = rhs[i];
  }
  return result;
}

Matrix hcat(const Matrix& lhs, const Vector& rhs)
{
  if (lhs.rows() != rhs.length()) {
    throw std::invalid_argument("Matrix rows must match vector length");
  }
  Matrix result(lhs.rows(), lhs.cols() + 1);
  for (size_t i = 0; i < lhs.rows(); ++i) {
    for (size_t j = 0; j < lhs.cols(); ++j) {
      result[i * result.cols() + j] = lhs[i * lhs.cols() + j];
    }
    result[i * result.cols() + lhs.cols()] = rhs[i];
  }
  return result;
}

Matrix hcat(const Vector& lhs, const Matrix& rhs)
{
  if (lhs.length() != rhs.rows()) {
    throw std::invalid_argument("Vector length must match matrix rows");
  }
  Matrix result(rhs.rows(), 1 + rhs.cols());
  for (size_t i = 0; i < rhs.rows(); ++i) {
    result[i * result.cols()] = lhs[i];
    for (size_t j = 0; j < rhs.cols(); ++j) {
      result[i * result.cols() + 1 + j] = rhs[i * rhs.cols() + j];
    }
  }
  return result;
}

Matrix vcat(const Matrix& lhs, const Matrix& rhs)
{
  if (lhs.cols() != rhs.cols()) {
    throw std::invalid_argument("Matrices must have the same number of columns");
  }
  Matrix result(lhs.rows() + rhs.rows(), lhs.cols());
  for (size_t i = 0; i < lhs.length(); ++i) {
    result[i] = lhs[i];
  }
  for (size_t i = 0; i < rhs.length(); ++i) {
    result[lhs.length() + i] = rhs[i];
  }
  return result;
}

real_t norm(const Matrix& mat, real_t p)
{
  if (p == 2.0) {
    real_t sum = 0.0;
    for (size_t i = 0; i < mat.length(); ++i) {
      sum += mat[i] * mat[i];
    }
    return std::sqrt(sum);
  }
  // For other p-norms, compute element-wise
  real_t sum = 0.0;
  for (size_t i = 0; i < mat.length(); ++i) {
    sum += std::pow(std::abs(mat[i]), p);
  }
  return std::pow(sum, 1.0 / p);
}

Matrix outer(const Vector& lhs, const Vector& rhs)
{
  Matrix result(lhs.length(), rhs.length());
  for (size_t i = 0; i < lhs.length(); ++i) {
    for (size_t j = 0; j < rhs.length(); ++j) {
      result[i * rhs.length() + j] = lhs[i] * rhs[j];
    }
  }
  return result;
}

Matrix inv(const Matrix& mat)
{
  if (mat.rows() != mat.cols()) {
    throw std::invalid_argument("Matrix must be square");
  }

  size_t n = mat.rows();
  Matrix augmented(n, 2 * n);

  // Create augmented matrix [A | I]
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      augmented[i * 2 * n + j] = mat[i * n + j];
    }
    augmented[i * 2 * n + n + i] = 1.0;
  }

  // Gaussian elimination with partial pivoting
  for (size_t k = 0; k < n; ++k) {
    // Find pivot
    size_t pivot_row = k;
    real_t max_val = std::abs(augmented[k * 2 * n + k]);
    for (size_t i = k + 1; i < n; ++i) {
      real_t val = std::abs(augmented[i * 2 * n + k]);
      if (val > max_val) {
        max_val = val;
        pivot_row = i;
      }
    }

    if (max_val < 1e-10) {
      throw std::invalid_argument("Matrix is singular");
    }

    // Swap rows
    if (pivot_row != k) {
      for (size_t j = 0; j < 2 * n; ++j) {
        std::swap(augmented[k * 2 * n + j], augmented[pivot_row * 2 * n + j]);
      }
    }

    // Scale pivot row
    real_t pivot = augmented[k * 2 * n + k];
    for (size_t j = 0; j < 2 * n; ++j) {
      augmented[k * 2 * n + j] /= pivot;
    }

    // Eliminate column
    for (size_t i = 0; i < n; ++i) {
      if (i != k) {
        real_t factor = augmented[i * 2 * n + k];
        for (size_t j = 0; j < 2 * n; ++j) {
          augmented[i * 2 * n + j] -= factor * augmented[k * 2 * n + j];
        }
      }
    }
  }

  // Extract inverse from augmented matrix
  Matrix result(n, n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      result[i * n + j] = augmented[i * 2 * n + n + j];
    }
  }

  return result;
}

}  // namespace sim