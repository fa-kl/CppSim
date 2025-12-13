/*****************************************************************************************
 * @file: Transform.hpp
 *
 * @brief: This file provides 3D rigid body transformations.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include <array>

#include "Matrix.hpp"
#include "Rotation.hpp"
#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

/**
 * @brief A class for handling 3D rigid body transformations (rotation + translation).
 *
 * @details Represents a rigid body transformation combining a rotation and a translation.
 * The transformation is applied as: T(v) = R * v + t, where R is the rotation and t is
 * the translation vector.
 */
class Transform
{
protected:
  /**
   * @brief The rotation component of the transformation.
   */
  Rotation m_rotation;

  /**
   * @brief The translation component of the transformation.
   */
  Vector m_translation;

public:
  /**
   * @brief Default constructor: identity transformation.
   */
  Transform();

  /**
   * @brief Construct a transformation from rotation and translation.
   * @param R Rotation component.
   * @param t Translation vector (must be 3-dimensional).
   * @throws std::invalid_argument If translation vector is not 3-dimensional.
   */
  Transform(const Rotation& R, const Vector& t);

  /**
   * @brief Get the 4x4 homogeneous transformation matrix.
   * @return 4x4 matrix representing the transformation.
   */
  Matrix getMatrix(void) const;

  /**
   * @brief Get the Euler angles of the rotation component.
   * @return The Euler angles (roll, pitch, yaw).
   */
  EulerAngles eulerAngles() const;

  /**
   * @brief Get the translation component.
   * @return The translation vector.
   */
  Vector translation() const;

  /**
   * @brief Get the rotation component.
   * @return The rotation.
   */
  Rotation rotation() const;

  /**
   * @brief Get the unit vectors of the rotated frame.
   * @return Three 3-dimensional unit vectors.
   */
  std::array<Vector, 3UL> unitVectors() const;

  /**
   * @brief Apply transformation to a vector.
   * @param T Transformation.
   * @param v 3-dimensional vector.
   * @return The transformed vector: R * v + t.
   */
  friend Vector operator*(const Transform& T, const Vector& v);

  /**
   * @brief Compose two transformations.
   * @param T1 First transformation.
   * @param T2 Second transformation.
   * @return The composed transformation T1 * T2 (apply T2 first, then T1).
   */
  friend Transform operator*(const Transform& T1, const Transform& T2);

  /**
   * @brief Apply a rotation after a transformation.
   * @param T Transformation.
   * @param R Rotation.
   * @return The composed transformation T * R.
   */
  friend Transform operator*(const Transform& T, const Rotation& R);

  /**
   * @brief Apply a transformation after a rotation.
   * @param R Rotation.
   * @param T Transformation.
   * @return The composed transformation R * T.
   */
  friend Transform operator*(const Rotation& R, const Transform& T);

  /**
   * @brief Compute the inverse transformation.
   * @param T Transformation.
   * @return The inverse transformation.
   */
  friend Transform inv(const Transform& T);

  /**
   * @brief Output transformation matrix to stream.
   * @param os Output stream.
   * @param T Transformation.
   * @return Reference to the output stream.
   */
  friend std::ostream& operator<<(std::ostream& os, const Transform& T);
};

/**
 * @brief Get the Euler angles from a transformation.
 * @param T Transformation.
 * @return The Euler angles (roll, pitch, yaw).
 */
inline EulerAngles eulerAngles(const Transform& T)
{
  return T.eulerAngles();
}

/**
 * @brief Get the translation vector from a transformation.
 * @param T Transformation.
 * @return The translation vector.
 */
inline Vector translation(const Transform& T)
{
  return T.translation();
}

/**
 * @brief Get the rotation from a transformation.
 * @param T Transformation.
 * @return The rotation component.
 */
inline Rotation rotation(const Transform& T)
{
  return T.rotation();
}

/**
 * @brief Get the unit vectors of the rotated frame.
 * @param T Transformation.
 * @return Three 3-dimensional unit vectors.
 */
inline std::array<Vector, 3UL> unitVectors(const Transform& T)
{
  return T.unitVectors();
}

}  // namespace sim
