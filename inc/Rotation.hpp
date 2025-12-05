/*****************************************************************************************
 * @file: Rotation.hpp
 *
 * @brief: This file provides rotation matrices in the form of a templated Rotation class.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "Matrix.hpp"
#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

template <size_t N>
class Transform;

template <size_t N>
class Rotation
{
protected:
  Matrix m_matrix;

public:
  Rotation() = default;

  // 2D constructor: angle
  template <size_t M = N>
  explicit Rotation(real_t angle, typename std::enable_if<M == 2>::type* = 0)
      : m_matrix({{std::cos(angle), -std::sin(angle)}, {std::sin(angle), std::cos(angle)}})
  {
  }

  // 3D constructor: roll, pitch, yaw
  template <size_t M = N>
  Rotation(real_t roll, real_t pitch, real_t yaw, typename std::enable_if<M == 3>::type* = 0)
  {
    m_matrix = (Rz(yaw) * Ry(pitch) * Rx(roll)).m_matrix;
  }

  explicit Rotation(const Matrix& mat)
  {
    dimension_t expected = {N, N};
    if (size(mat) != expected) {
      throw std::invalid_argument("Incompatible dimension for rotation matrix");
    }
    m_matrix = mat;
  }

  // 2D: single angle
  template <size_t M = N>
  typename std::enable_if<M == 2, real_t>::type angle() const
  {
    return std::atan2(m_matrix(2, 1), m_matrix(1, 1));
  }

  // 3D: return Euler angles as Vector (roll, pitch, yaw)
  template <size_t M = N>
  typename std::enable_if<M == 3, Vector>::type eulerAngles() const
  {
    real_t sy = -m_matrix(3, 1);
    real_t cy = std::sqrt(m_matrix(1, 1) * m_matrix(1, 1) + m_matrix(2, 1) * m_matrix(2, 1));
    real_t yaw = std::atan2(m_matrix(2, 1), m_matrix(1, 1));
    real_t pitch = std::atan2(sy, cy);
    real_t roll = std::atan2(m_matrix(3, 2), m_matrix(3, 3));
    return Vector({roll, pitch, yaw});
  }

  std::vector<Vector> unitVectors() const
  {
    std::vector<Vector> res;
    for (index_t c = 1; c <= static_cast<index_t>(N); ++c) {
      Vector col(N);
      for (index_t r = 1; r <= static_cast<index_t>(N); ++r) {
        col(r) = m_matrix(r, c);
      }
      res.push_back(col);
    }
    return res;
  }

  friend Vector operator*(const Rotation& R, const Vector& v) { return R.m_matrix * v; }

  friend Rotation operator*(const Rotation& R1, const Rotation& R2) { return Rotation(R1.m_matrix * R2.m_matrix); }

  friend Rotation inv(const Rotation& R) { return Rotation(transpose(R.m_matrix)); }

  template <size_t M = N>
  static typename std::enable_if<M == 3, Rotation>::type Rx(real_t roll)
  {
    Matrix m({{1, 0, 0}, {0, std::cos(roll), -std::sin(roll)}, {0, std::sin(roll), std::cos(roll)}});
    return Rotation(m);
  }

  template <size_t M = N>
  static typename std::enable_if<M == 3, Rotation>::type Ry(real_t pitch)
  {
    Matrix m({{std::cos(pitch), 0, std::sin(pitch)}, {0, 1, 0}, {-std::sin(pitch), 0, std::cos(pitch)}});
    return Rotation(m);
  }

  template <size_t M = N>
  static typename std::enable_if<M == 3, Rotation>::type Rz(real_t yaw)
  {
    Matrix m({{std::cos(yaw), -std::sin(yaw), 0}, {std::sin(yaw), std::cos(yaw), 0}, {0, 0, 1}});
    return Rotation(m);
  }
};

inline real_t angle(const Rotation<2>& R)
{
  return R.angle();
}

inline Vector eulerAngles(const Rotation<3>& R)
{
  return R.eulerAngles();
}

inline std::vector<Vector> unitVectors(const Rotation<2>& R)
{
  return R.unitVectors();
}

inline std::vector<Vector> unitVectors(const Rotation<3>& R)
{
  return R.unitVectors();
}

}  // namespace sim
