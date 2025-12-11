/*****************************************************************************************
 * @file: Vector.cpp
 *
 * @brief: This file implements the Vector's functions.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Vector.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Matrix.hpp"

namespace sim
{

/*****************************************************************************************
 * Constructors and Destructors
 ****************************************************************************************/

Vector::Vector() : m_data() {}

Vector::Vector(size_t n) : m_data(n, 0.0) {}

Vector::Vector(std::initializer_list<real_t> values) : m_data(values) {}

Vector::Vector(const Vector& other) : m_data(other.m_data) {}

Vector::Vector(Vector&& other) noexcept : m_data(std::move(other.m_data)) {}

Vector::~Vector() = default;

/*****************************************************************************************
 * Assignment Operators
 ****************************************************************************************/

Vector& Vector::operator=(const Vector& other)
{
  if (this != &other) {
    m_data = other.m_data;
  }
  return *this;
}

Vector& Vector::operator=(Vector&& other) noexcept
{
  if (this != &other) {
    m_data = std::move(other.m_data);
  }
  return *this;
}

/*****************************************************************************************
 * Size and Capacity
 ****************************************************************************************/

size_t Vector::length() const
{
  return m_data.size();
}

dimension_t Vector::size() const
{
  return {m_data.size(), 1};
}

size_t Vector::rows() const
{
  return m_data.size();
}

size_t Vector::cols() const
{
  return 1;
}

bool_t Vector::isEmpty() const
{
  return m_data.empty();
}

/***************************************************************************************
 * Element Access
 **************************************************************************************/

real_t& Vector::operator[](size_t i)
{
  return m_data[i];
}

const real_t& Vector::operator[](size_t i) const
{
  return m_data[i];
}

real_t& Vector::operator()(index_t i)
{
  if (i == 0) {
    throw std::out_of_range("Index cannot be zero (1-based indexing)");
  }
  size_t idx = (i > 0) ? static_cast<size_t>(i - 1) : (m_data.size() + static_cast<size_t>(i));
  if (idx >= m_data.size()) {
    throw std::out_of_range("Index out of range");
  }
  return m_data[idx];
}

const real_t& Vector::operator()(index_t i) const
{
  if (i == 0) {
    throw std::out_of_range("Index cannot be zero (1-based indexing)");
  }
  size_t idx = (i > 0) ? static_cast<size_t>(i - 1) : (m_data.size() + static_cast<size_t>(i));
  if (idx >= m_data.size()) {
    throw std::out_of_range("Index out of range");
  }
  return m_data[idx];
}

/*****************************************************************************************
 * Unary Operators
 ****************************************************************************************/

Vector Vector::operator+() const
{
  return *this;
}

Vector Vector::operator-() const
{
  Vector result(m_data.size());
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = -m_data[i];
  }
  return result;
}

/*****************************************************************************************
 * Arithmetic Operators
 ****************************************************************************************/

Vector Vector::operator+(const Vector& rhs) const
{
  if (m_data.size() != rhs.m_data.size()) {
    throw std::invalid_argument("Incompatible vector sizes");
  }
  Vector result(m_data.size());
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = m_data[i] + rhs.m_data[i];
  }
  return result;
}

Vector Vector::operator+(real_t scalar) const
{
  Vector result(m_data.size());
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = m_data[i] + scalar;
  }
  return result;
}

Vector Vector::operator-(const Vector& rhs) const
{
  if (m_data.size() != rhs.m_data.size()) {
    throw std::invalid_argument("Incompatible vector sizes");
  }
  Vector result(m_data.size());
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = m_data[i] - rhs.m_data[i];
  }
  return result;
}

Vector Vector::operator-(real_t scalar) const
{
  Vector result(m_data.size());
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = m_data[i] - scalar;
  }
  return result;
}

real_t Vector::operator*(const Vector& rhs) const
{
  if (m_data.size() != rhs.m_data.size()) {
    throw std::invalid_argument("Incompatible vector sizes");
  }
  real_t result = 0.0;
  for (size_t i = 0; i < m_data.size(); ++i) {
    result += m_data[i] * rhs.m_data[i];
  }
  return result;
}

Vector Vector::operator*(real_t scalar) const
{
  Vector result(m_data.size());
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = m_data[i] * scalar;
  }
  return result;
}

Vector Vector::operator/(real_t scalar) const
{
  if (std::abs(scalar) < 1e-10) {
    throw std::invalid_argument("Division by zero");
  }
  Vector result(m_data.size());
  for (size_t i = 0; i < m_data.size(); ++i) {
    result.m_data[i] = m_data[i] / scalar;
  }
  return result;
}

/**
 * @brief Multiply a scalar by a vector.
 *
 * @param scalar Scalar value.
 * @param vec Vector.
 * @return Result vector.
 */
Vector operator*(real_t scalar, const Vector& vec)
{
  return vec * scalar;
}

/**
 * @brief Add a scalar to a vector.
 *
 * @param scalar Scalar value.
 * @param vec Vector.
 * @return Result vector.
 */
Vector operator+(real_t scalar, const Vector& vec)
{
  return vec + scalar;
}

/**
 * @brief Subtract a vector from a scalar.
 *
 * @param scalar Scalar value.
 * @param vec Vector.
 * @return Result vector.
 */
Vector operator-(real_t scalar, const Vector& vec)
{
  Vector result(vec.m_data.size());
  for (size_t i = 0; i < vec.m_data.size(); ++i) {
    result.m_data[i] = scalar - vec.m_data[i];
  }
  return result;
}

/**
 * @brief Divide a scalar by all elements of a vector.
 *
 * @param scalar Scalar value (dividend).
 * @param vec Vector (divisor).
 * @return Result vector.
 */
Vector operator/(real_t scalar, const Vector& vec)
{
  Vector result(vec.m_data.size());
  for (size_t i = 0; i < vec.m_data.size(); ++i) {
    if (std::abs(vec.m_data[i]) < 1e-10) {
      throw std::invalid_argument("Division by zero");
    }
    result.m_data[i] = scalar / vec.m_data[i];
  }
  return result;
}

/*****************************************************************************************
 * Comparison Operators
 ****************************************************************************************/

bool_t Vector::operator==(const Vector& rhs) const
{
  if (m_data.size() != rhs.m_data.size()) {
    return false;
  }
  for (size_t i = 0; i < m_data.size(); ++i) {
    if (m_data[i] != rhs.m_data[i]) {
      return false;
    }
  }
  return true;
}

bool_t Vector::operator!=(const Vector& rhs) const
{
  return !(*this == rhs);
}

/*****************************************************************************************
 * Print Function(s)
 ****************************************************************************************/

/**
 * @brief Output stream operator for vectors.
 *
 * @param os Output stream.
 * @param vec Vector to output.
 * @return Output stream.
 */
std::ostream& operator<<(std::ostream& os, const Vector& vec)
{
  // Format each element with fixed 3 decimal places and align columns
  std::vector<std::string> parts;
  parts.reserve(vec.length());

  std::ostringstream tmp;
  tmp.setf(std::ios::fixed, std::ios::floatfield);
  tmp.precision(3);

  size_t max_width = 0;
  for (size_t i = 0; i < vec.length(); ++i) {
    tmp.str("");
    tmp.clear();
    tmp << vec[i];
    std::string s = tmp.str();
    parts.push_back(s);
    if (s.size() > max_width)
      max_width = s.size();
  }

  os << "[";
  for (size_t i = 0; i < parts.size(); ++i) {
    os << std::setw(static_cast<int>(max_width)) << parts[i];
    if (i < parts.size() - 1) {
      os << "\n ";
    }
  }
  os << "]";
  return os;
}

/*****************************************************************************************
 * Non-Member Vector Functions
 ****************************************************************************************/

Vector zeros(size_t n)
{
  return Vector(n);
}

Vector ones(size_t n)
{
  Vector result(n);
  for (size_t i = 0; i < n; ++i) {
    result[i] = 1.0;
  }
  return result;
}

size_t length(const Vector& vec)
{
  return vec.length();
}

dimension_t size(const Vector& vec)
{
  return vec.size();
}

bool_t isEmpty(const Vector& vec)
{
  return vec.isEmpty();
}

real_t norm(const Vector& vec, real_t p)
{
  if (vec.isEmpty()) {
    return 0.0;
  }

  // Infinity norm
  if (std::isinf(p) && p > 0) {
    real_t max_val = 0.0;
    for (size_t i = 0; i < vec.length(); ++i) {
      max_val = std::max(max_val, std::abs(vec[i]));
    }
    return max_val;
  }

  // -Infinity norm
  if (std::isinf(p) && p < 0) {
    real_t min_val = std::numeric_limits<real_t>::infinity();
    for (size_t i = 0; i < vec.length(); ++i) {
      min_val = std::min(min_val, std::abs(vec[i]));
    }
    return min_val;
  }

  if (p <= 0) {
    throw std::invalid_argument("Norm parameter p must be positive");
  }

  // 1-norm
  if (p == 1.0) {
    real_t sum = 0.0;
    for (size_t i = 0; i < vec.length(); ++i) {
      sum += std::abs(vec[i]);
    }
    return sum;
  }

  // 2-norm (Euclidean)
  if (p == 2.0) {
    real_t sum = 0.0;
    for (size_t i = 0; i < vec.length(); ++i) {
      sum += vec[i] * vec[i];
    }
    return std::sqrt(sum);
  }

  // p-norm
  real_t sum = 0.0;
  for (size_t i = 0; i < vec.length(); ++i) {
    sum += std::pow(std::abs(vec[i]), p);
  }
  return std::pow(sum, 1.0 / p);
}

Vector normalize(const Vector& vec)
{
  real_t n = norm(vec);
  if (std::abs(n) < 1e-10) {
    throw std::invalid_argument("Cannot normalize zero vector");
  }
  return vec / n;
}

Vector abs(const Vector& vec)
{
  Vector result(vec.length());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = std::abs(vec[i]);
  }
  return result;
}

Vector exp(const Vector& vec)
{
  Vector result(vec.length());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = std::exp(vec[i]);
  }
  return result;
}

Vector log(const Vector& vec)
{
  Vector result(vec.length());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = std::log(vec[i]);
  }
  return result;
}

Vector sin(const Vector& vec)
{
  Vector result(vec.length());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = std::sin(vec[i]);
  }
  return result;
}

Vector cos(const Vector& vec)
{
  Vector result(vec.length());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = std::cos(vec[i]);
  }
  return result;
}

Vector tan(const Vector& vec)
{
  Vector result(vec.length());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = std::tan(vec[i]);
  }
  return result;
}

Vector asin(const Vector& vec)
{
  Vector result(vec.length());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = std::asin(vec[i]);
  }
  return result;
}

Vector acos(const Vector& vec)
{
  Vector result(vec.length());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = std::acos(vec[i]);
  }
  return result;
}

Vector atan(const Vector& vec)
{
  Vector result(vec.length());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = std::atan(vec[i]);
  }
  return result;
}

Vector pow(const Vector& vec, real_t exponent)
{
  Vector result(vec.length());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = std::pow(vec[i], exponent);
  }
  return result;
}

Vector sqrt(const Vector& vec)
{
  Vector result(vec.length());
  for (size_t i = 0; i < vec.length(); ++i) {
    result[i] = std::sqrt(vec[i]);
  }
  return result;
}

real_t sum(const Vector& vec)
{
  real_t result = 0.0;
  for (size_t i = 0; i < vec.length(); ++i) {
    result += vec[i];
  }
  return result;
}

Vector diff(const Vector& vec)
{
  if (vec.length() < 2) {
    return Vector(0);
  }
  Vector result(vec.length() - 1);
  for (size_t i = 0; i < vec.length() - 1; ++i) {
    result[i] = vec[i + 1] - vec[i];
  }
  return result;
}

Vector elmul(const Vector& lhs, const Vector& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw std::invalid_argument("Incompatible vector sizes");
  }
  Vector result(lhs.length());
  for (size_t i = 0; i < lhs.length(); ++i) {
    result[i] = lhs[i] * rhs[i];
  }
  return result;
}

Vector eldiv(const Vector& lhs, const Vector& rhs)
{
  if (lhs.length() != rhs.length()) {
    throw std::invalid_argument("Incompatible vector sizes");
  }
  Vector result(lhs.length());
  for (size_t i = 0; i < lhs.length(); ++i) {
    if (std::abs(rhs[i]) < 1e-10) {
      throw std::invalid_argument("Division by zero");
    }
    result[i] = lhs[i] / rhs[i];
  }
  return result;
}

real_t dot(const Vector& lhs, const Vector& rhs)
{
  return lhs * rhs;
}

Vector cross(const Vector& lhs, const Vector& rhs)
{
  if (lhs.length() != 3 || rhs.length() != 3) {
    throw std::invalid_argument("Cross product requires 3D vectors");
  }
  Vector result(3);
  result[0] = lhs[1] * rhs[2] - lhs[2] * rhs[1];
  result[1] = lhs[2] * rhs[0] - lhs[0] * rhs[2];
  result[2] = lhs[0] * rhs[1] - lhs[1] * rhs[0];
  return result;
}

real_t max(const Vector& vec)
{
  if (vec.isEmpty()) {
    throw std::invalid_argument("Cannot find max of empty vector");
  }
  real_t max_val = vec[0];
  for (size_t i = 1; i < vec.length(); ++i) {
    max_val = std::max(max_val, vec[i]);
  }
  return max_val;
}

real_t min(const Vector& vec)
{
  if (vec.isEmpty()) {
    throw std::invalid_argument("Cannot find min of empty vector");
  }
  real_t min_val = vec[0];
  for (size_t i = 1; i < vec.length(); ++i) {
    min_val = std::min(min_val, vec[i]);
  }
  return min_val;
}

real_t mean(const Vector& vec)
{
  if (vec.isEmpty()) {
    return 0.0;
  }
  return sum(vec) / static_cast<real_t>(vec.length());
}

real_t var(const Vector& vec)
{
  if (vec.isEmpty()) {
    return 0.0;
  }
  real_t m = mean(vec);
  real_t sum_sq = 0.0;
  for (size_t i = 0; i < vec.length(); ++i) {
    real_t diff = vec[i] - m;
    sum_sq += diff * diff;
  }
  return sum_sq / static_cast<real_t>(vec.length());
}

real_t std(const Vector& vec)
{
  return std::sqrt(var(vec));
}

Vector vcat(const Vector& lhs, const Vector& rhs)
{
  Vector result(lhs.length() + rhs.length());
  for (size_t i = 0; i < lhs.length(); ++i) {
    result[i] = lhs[i];
  }
  for (size_t i = 0; i < rhs.length(); ++i) {
    result[lhs.length() + i] = rhs[i];
  }
  return result;
}

}  // namespace sim