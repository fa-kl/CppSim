/*****************************************************************************************
 * @file: Rotation.cpp
 *
 * @brief: This file implements the Rotation class functions.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Rotation.hpp"

#include <cmath>
#include <stdexcept>

#include "Matrix.hpp"
#include "Transform.hpp"
#include "Vector.hpp"

namespace sim
{

/*****************************************************************************************
 * Constructors
 ****************************************************************************************/

Rotation::Rotation(real_t angle) : m_matrix({{std::cos(angle), -std::sin(angle)}, {std::sin(angle), std::cos(angle)}})
{
}

Rotation::Rotation(const Matrix& mat)
{
  if (size(mat) != dimension_t{2, 2}) {
    throw std::invalid_argument("Incompatible dimension for rotation matrix");
  }
  m_matrix = mat;
}

/*****************************************************************************************
 * Member Functions
 ****************************************************************************************/

real_t Rotation::angle() const
{
  return std::atan2(m_matrix(2, 1), m_matrix(1, 1));
}

Matrix Rotation::unitVectors() const
{
  return m_matrix;
}

/*****************************************************************************************
 * Non-Member Functions
 ****************************************************************************************/

real_t angle(const Rotation& R)
{
  return R.angle();
}

Matrix unitVectors(const Rotation& R)
{
  return R.unitVectors();
}

Vector operator*(const Rotation& R, const Vector& v)
{
  return R.m_matrix * v;
}

Rotation operator*(const Rotation& R1, const Rotation& R2)
{
  return Rotation(R1.m_matrix * R2.m_matrix);
}

Rotation inv(const Rotation& R)
{
  return Rotation(transpose(R.m_matrix));
}

}  // namespace sim
