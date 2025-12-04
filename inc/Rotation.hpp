/*****************************************************************************************
 * @file: Rotation.hpp
 *
 * @brief: This file provides rotation matrices in the form of a Rotation class.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include "Matrix.hpp"
#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

class Transform;

/**
 * @brief A class representing 2D rotations.
 *
 * @details This class encapsulates 2D rotation transformations using a 2x2 rotation matrix.
 * Rotations can be composed, inverted, and applied to vectors. The rotation angle is measured
 * in radians, counter-clockwise from the positive x-axis.
 */
class Rotation
{
protected:
  /**
   * @brief The 2x2 rotation matrix.
   */
  Matrix m_matrix;

public:
  /**
   * @brief Construct a rotation from an angle.
   *
   * @param angle Rotation angle in radians (counter-clockwise).
   */
  explicit Rotation(real_t angle);

  /**
   * @brief Construct a rotation from a 2x2 matrix.
   *
   * @param mat 2x2 rotation matrix.
   * @throws std::invalid_argument If matrix is not 2x2.
   */
  explicit Rotation(const Matrix& mat);

  /**
   * @brief Get the rotation angle.
   *
   * @return Rotation angle in radians (counter-clockwise).
   */
  real_t angle() const;

  /**
   * @brief Get the unit vectors defining this rotation.
   *
   * @details Returns the rotated standard basis vectors [e1, e2] where e1 is the
   * rotated x-axis and e2 is the rotated y-axis.
   *
   * @return Matrix with columns representing the rotated unit vectors.
   */
  Matrix unitVectors() const;

  friend Vector operator*(const Rotation& R, const Vector& v);
  friend Rotation operator*(const Rotation& R1, const Rotation& R2);
  friend Transform operator*(const Transform& T, const Rotation& R);
  friend Transform operator*(const Rotation& R, const Transform& T);
  friend Rotation inv(const Rotation& R);
  friend Matrix unitVectors(const Rotation& R);
};

/**
 * @brief Get the rotation angle.
 *
 * @param R Rotation object.
 * @return Rotation angle in radians (counter-clockwise).
 */
real_t angle(const Rotation& R);

/**
 * @brief Get the unit vectors defining the rotation.
 *
 * @details Returns the rotated standard basis vectors [e1, e2] where e1 is the
 * rotated x-axis and e2 is the rotated y-axis.
 *
 * @param R Rotation object.
 * @return Matrix with columns representing the rotated unit vectors.
 */
Matrix unitVectors(const Rotation& R);

/**
 * @brief Apply rotation to a vector.
 *
 * @param R Rotation object.
 * @param v Vector to rotate.
 * @return Rotated vector.
 */
Vector operator*(const Rotation& R, const Vector& v);

/**
 * @brief Compose two rotations.
 *
 * @details Computes R1 * R2, which applies R2 first, then R1.
 *
 * @param R1 First rotation.
 * @param R2 Second rotation.
 * @return Composed rotation.
 */
Rotation operator*(const Rotation& R1, const Rotation& R2);

/**
 * @brief Compute the inverse of a rotation.
 *
 * @details The inverse of a rotation matrix is its transpose.
 *
 * @param R Rotation object.
 * @return Inverse rotation.
 */
Rotation inv(const Rotation& R);

}  // namespace sim