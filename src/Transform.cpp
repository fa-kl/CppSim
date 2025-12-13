/*****************************************************************************************
 * @file: Transform.cpp
 *
 * @brief: Implementation of 3D rigid body transformations.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Transform.hpp"

#include <stdexcept>

namespace sim
{

Transform::Transform() : m_rotation(Rotation()), m_translation(zeros(3)) {}

Transform::Transform(const Rotation& R, const Vector& t) : m_rotation(R)
{
  if (length(t) != 3) {
    throw std::invalid_argument("Translation vector must match transform dimension");
  }
  m_translation = t;
}

Matrix Transform::getMatrix(void) const
{
  return vcat(hcat(m_rotation.getMatrix(), m_translation), {{0, 0, 0, 1}});
}

EulerAngles Transform::eulerAngles() const
{
  return m_rotation.eulerAngles();
}

Vector Transform::translation() const
{
  return m_translation;
}

Rotation Transform::rotation() const
{
  return m_rotation;
}

std::array<Vector, 3UL> Transform::unitVectors() const
{
  return m_rotation.unitVectors();
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

std::ostream& operator<<(std::ostream& os, const Transform& T)
{
  os << T.getMatrix();
  return os;
}

}  // namespace sim
