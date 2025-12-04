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

#include "Rotation.hpp"
#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

/**
 * @brief A class representing 2D rigid body transformations (rotation + translation).
 *
 * @details This class represents a 2D rigid body transformation consisting of a rotation
 * followed by a translation. Transforms can be composed, inverted, and applied to vectors.
 * The transformation T(v) = R*v + t where R is rotation and t is translation.
 */
class Transform
{
protected:
  /**
   * @brief The rotation component.
   */
  Rotation m_rotation;

  /**
   * @brief The translation component.
   */
  Vector m_translation;

public:
  /**
   * @brief Construct a transform from rotation and translation.
   *
   * @param R Rotation component.
   * @param t Translation vector (must be 2-dimensional).
   * @throws std::invalid_argument If translation vector is not 2-dimensional.
   */
  Transform(const Rotation& R, const Vector& t);

  /**
   * @brief Get the rotation angle.
   *
   * @return Rotation angle in radians (counter-clockwise).
   */
  real_t angle() const;

  /**
   * @brief Get the translation vector.
   *
   * @return Translation vector.
   */
  Vector translation() const;

  /**
   * @brief Get the rotation component.
   *
   * @return Rotation object.
   */
  Rotation rotation() const;

  /**
   * @brief Get the unit vectors defining the rotation component.
   *
   * @details Returns the rotated standard basis vectors [e1, e2] where e1 is the
   * rotated x-axis and e2 is the rotated y-axis.
   *
   * @return Matrix with columns representing the rotated unit vectors.
   */
  Matrix unitVectors() const;

  friend Vector operator*(const Transform& T, const Vector& v);
  friend Transform operator*(const Transform& T1, const Transform& T2);
  friend Transform operator*(const Transform& T, const Rotation& R);
  friend Transform operator*(const Rotation& R, const Transform& T);
  friend Transform inv(const Transform& T);
  friend Matrix unitVectors(const Transform& T);
};

/**
 * @brief Get the rotation angle of a transform.
 *
 * @param T Transform object.
 * @return Rotation angle in radians (counter-clockwise).
 */
real_t angle(const Transform& T);

/**
 * @brief Get the translation vector of a transform.
 *
 * @param T Transform object.
 * @return Translation vector.
 */
Vector translation(const Transform& T);

/**
 * @brief Get the rotation component of a transform.
 *
 * @param T Transform object.
 * @return Rotation object.
 */
Rotation rotation(const Transform& T);

/**
 * @brief Get the unit vectors defining the rotation component of a transform.
 *
 * @details Returns the rotated standard basis vectors [e1, e2] where e1 is the
 * rotated x-axis and e2 is the rotated y-axis.
 *
 * @param T Transform object.
 * @return Matrix with columns representing the rotated unit vectors.
 */
Matrix unitVectors(const Transform& T);

/**
 * @brief Apply transform to a vector.
 *
 * @details Computes T(v) = R*v + t where R is rotation and t is translation.
 *
 * @param T Transform object.
 * @param v Vector to transform.
 * @return Transformed vector.
 */
Vector operator*(const Transform& T, const Vector& v);

/**
 * @brief Compose two transforms.
 *
 * @details Computes T1 * T2, which applies T2 first, then T1.
 * Mathematically: (T1 * T2)(v) = T1(T2(v)) = R1*(R2*v + t2) + t1 = (R1*R2)*v + (R1*t2 + t1).
 *
 * @param T1 First transform.
 * @param T2 Second transform.
 * @return Composed transform.
 */
Transform operator*(const Transform& T1, const Transform& T2);

/**
 * @brief Compose a transform with a rotation.
 *
 * @details Computes T * R, which applies R first, then T.
 *
 * @param T Transform object.
 * @param R Rotation object.
 * @return Composed transform.
 */
Transform operator*(const Transform& T, const Rotation& R);

/**
 * @brief Compose a rotation with a transform.
 *
 * @details Computes R * T, which applies T first, then R.
 *
 * @param R Rotation object.
 * @param T Transform object.
 * @return Composed transform.
 */
Transform operator*(const Rotation& R, const Transform& T);

/**
 * @brief Compute the inverse of a transform.
 *
 * @details The inverse transform T^(-1) satisfies T^(-1)(T(v)) = v for all v.
 * If T(v) = R*v + t, then T^(-1)(v) = R^(-1)*(v - t) = R^(-1)*v - R^(-1)*t.
 *
 * @param T Transform object.
 * @return Inverse transform.
 */
Transform inv(const Transform& T);

}  // namespace sim
