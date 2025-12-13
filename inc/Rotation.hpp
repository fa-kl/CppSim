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

#include <array>

#include "Matrix.hpp"
#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

/**
 * @brief A struct type for euler angles of 3D rotations.
 *
 * @param roll Angle of rotation around the x-axis.
 * @param pitch Angle of rotation around the y-axis.
 * @param yaw Angle of rotation around the z-axis.
 */
struct EulerAngles
{
    real_t roll;
    real_t pitch;
    real_t yaw;
};

/**
 * @brief Convert radiants to degrees.
 */
inline real_t rad2deg(real_t rad)
{
    return rad * 180.0 / M_PI;
}

/**
 * @brief Convert degrees to radiants.
 */
inline real_t deg2rad(real_t deg)
{
    return deg / 180.0 * M_PI;
}

class Transform;

/**
 * @brief A class for handling 3D rotations.
 */
class Rotation
{
  protected:
    /**
     * @brief The internal rotation matrix.
     */
    Matrix m_matrix;

  public:
    /**
     * @brief Default constructor: unit rotation.
     */
    Rotation();

    /**
     * @brief Creates a rotation from roll, pitch, and yaw angles.
     * @param roll Angle of rotation around the x-axis in radiants.
     * @param pitch Angle of rotation around the y-axis in radiants.
     * @param yaw Angle of rotation around the z-axis in radiants.
     */
    Rotation(real_t roll, real_t pitch, real_t yaw);

    /**
     * @brief Creates a rotation from euler angles.
     * @param angles Euler angles.
     */
    Rotation(EulerAngles angles);

    /**
     * @brief Creates a rotation from a matrix.
     * @param mat 3x3 matrix.
     */
    explicit Rotation(const Matrix& mat);

    /**
     * @brief Get the underlying rotation matrix.
     * @return The 3x3 rotation matrix.
     */
    Matrix getMatrix() const;

    /**
     * @brief Get the Euler angles corresponding to the rotation.
     * @return The Euler angles.
     */
    EulerAngles eulerAngles() const;

    /**
     * @brief Get the unit vectors of the rotated frame.
     * @return Three 3-dimensional unit vectors.
     */
    std::array<Vector, 3UL> unitVectors() const;

    /**
     * @brief Rotates a vector.
     * @param R Rotation.
     * @param v 3-dimensional vector.
     * @return The rotated vector.
     */
    friend Vector operator*(const Rotation& R, const Vector& v);

    /**
     * @brief Compose two rotations.
     * @param R1 First rotation.
     * @param R2 Second rotation.
     * @return The composed rotation R1 * R2 (apply R2 first, then R1).
     */
    friend Rotation operator*(const Rotation& R1, const Rotation& R2);

    /**
     * @brief Compute the inverse rotation.
     * @param R Rotation.
     * @return The inverse rotation (transpose of rotation matrix).
     */
    friend Rotation inv(const Rotation& R);

    /**
     * @brief Create a rotation around the x-axis.
     * @param roll Angle of rotation in radians.
     * @return Rotation matrix for rotation around x-axis.
     */
    static Rotation Rx(real_t roll);

    /**
     * @brief Create a rotation around the y-axis.
     * @param pitch Angle of rotation in radians.
     * @return Rotation matrix for rotation around y-axis.
     */
    static Rotation Ry(real_t pitch);

    /**
     * @brief Create a rotation around the z-axis.
     * @param yaw Angle of rotation in radians.
     * @return Rotation matrix for rotation around z-axis.
     */
    static Rotation Rz(real_t yaw);

    /**
     * @brief Output rotation matrix to stream.
     * @param os Output stream.
     * @param R Rotation.
     * @return Reference to the output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const Rotation& R);
};

/**
 * @brief Get the Euler angles from a rotation.
 * @param R Rotation.
 * @return The Euler angles (roll, pitch, yaw).
 */
inline EulerAngles eulerAngles(const Rotation& R)
{
    return R.eulerAngles();
}

/**
 * @brief Get the unit vectors of the rotated frame.
 * @param R Rotation.
 * @return Three 3-dimensional unit vectors representing the rotated coordinate axes.
 */
inline std::array<Vector, 3UL> unitVectors(const Rotation& R)
{
    return R.unitVectors();
}

}  // namespace sim
