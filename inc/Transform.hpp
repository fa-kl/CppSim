/*****************************************************************************************
 * @file: Transform.hpp
 *
 * @brief: This file provides 2D rigid body transformations.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include <stdexcept>
#include <type_traits>
#include <vector>

#include "Matrix.hpp"
#include "Rotation.hpp"
#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

template <size_t N = 3>
class Transform
{
protected:
  Rotation<N> m_rotation;
  Vector m_translation;

public:
  Transform() : m_rotation(Rotation<N>()), m_translation(zeros(3)) {}

  Transform(const Rotation<N>& R, const Vector& t) : m_rotation(R)
  {
    if (length(t) != N) {
      throw std::invalid_argument("Translation vector must match transform dimension");
    }
    m_translation = t;
  }

  // 2D: scalar angle
  template <size_t M = N>
  typename std::enable_if<M == 2, real_t>::type angle() const
  {
    return m_rotation.angle();
  }

  // 3D: Euler angles
  template <size_t M = N>
  typename std::enable_if<M == 3, EulerAngles_t>::type eulerAngles() const
  {
    return m_rotation.eulerAngles();
  }

  Vector translation() const { return m_translation; }

  Rotation<N> rotation() const { return m_rotation; }

  std::vector<Vector> unitVectors() const { return m_rotation.unitVectors(); }

  friend Vector operator*(const Transform& T, const Vector& v) { return T.m_rotation * v + T.m_translation; }

  friend Transform operator*(const Transform& T1, const Transform& T2)
  {
    return Transform(T1.m_rotation * T2.m_rotation, T1.m_rotation * T2.m_translation + T1.m_translation);
  }

  friend Transform operator*(const Transform& T, const Rotation<N>& R)
  {
    return Transform(T.m_rotation * R, T.m_translation);
  }

  friend Transform operator*(const Rotation<N>& R, const Transform& T)
  {
    return Transform(R * T.m_rotation, R * T.m_translation);
  }

  friend Transform inv(const Transform& T)
  {
    Rotation<N> R_inv = inv(T.m_rotation);
    return Transform(R_inv, R_inv * (-T.m_translation));
  }

  friend std::ostream& operator<<(std::ostream& os, const Transform& T)
  {
    if constexpr (N == 2) {
      os << vcat(hcat(T.m_rotation.getMatrix(), T.m_translation), {{0, 0, 1}});
    } else if constexpr (N == 3) {
      os << vcat(hcat(T.m_rotation.getMatrix(), T.m_translation), {{0, 0, 0, 1}});
    } else {
      throw std::invalid_argument("Invalid N.");
    }
    return os;
  }
};

inline real_t angle(const Transform<2>& T)
{
  return T.angle();
}

inline EulerAngles_t eulerAngles(const Transform<3>& T)
{
  return T.eulerAngles();
}

inline Vector translation(const Transform<2>& T)
{
  return T.translation();
}

inline Vector translation(const Transform<3>& T)
{
  return T.translation();
}

inline Rotation<2> rotation(const Transform<2>& T)
{
  return T.rotation();
}

inline Rotation<3> rotation(const Transform<3>& T)
{
  return T.rotation();
}

inline std::vector<Vector> unitVectors(const Transform<2>& T)
{
  return T.unitVectors();
}

inline std::vector<Vector> unitVectors(const Transform<3>& T)
{
  return T.unitVectors();
}

}  // namespace sim
