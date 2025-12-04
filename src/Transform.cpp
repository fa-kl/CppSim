/*****************************************************************************************
 * @file: Transform.cpp
 *
 * @brief: This file implements the Transform class functions.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Transform.hpp"

#include <stdexcept>

#include "Rotation.hpp"
#include "Vector.hpp"

namespace sim
{

/*****************************************************************************************
 * Constructors
 ****************************************************************************************/

Transform::Transform(const Rotation& R, const Vector& t) : m_rotation(R)
{
  if (length(t) != 2) {
    throw std::invalid_argument("Translation vector must be 2-dimensional");
  }
  m_translation = t;
}

/*****************************************************************************************
 * Member Functions
 ****************************************************************************************/

real_t Transform::angle() const
{
  return m_rotation.angle();
}

Vector Transform::translation() const
{
  return m_translation;
}

Rotation Transform::rotation() const
{
  return m_rotation;
}

Matrix Transform::unitVectors() const
{
  return m_rotation.unitVectors();
}

/*****************************************************************************************
 * Non-Member Functions
 ****************************************************************************************/

real_t angle(const Transform& T)
{
  return T.angle();
}

Vector translation(const Transform& T)
{
  return T.translation();
}

Rotation rotation(const Transform& T)
{
  return T.rotation();
}

Matrix unitVectors(const Transform& T)
{
  return T.unitVectors();
}

Vector operator*(const Transform& T, const Vector& v)
{
  return T.m_rotation * v + T.m_translation;
}

Transform operator*(const Transform& T1, const Transform& T2)
{
  return Transform(T1.m_rotation * T2.m_rotation, T1.m_rotation * T2.m_translation + T1.m_translation);
}

Transform operator*(const Transform& T, const Rotation& R)
{
  return Transform(T.m_rotation * R, T.m_translation);
}

Transform operator*(const Rotation& R, const Transform& T)
{
  return Transform(R * T.m_rotation, R * T.m_translation);
}

Transform inv(const Transform& T)
{
  Rotation R_inv = inv(T.m_rotation);
  return Transform(R_inv, R_inv * (-T.m_translation));
}

}  // namespace sim
