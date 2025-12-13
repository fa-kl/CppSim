/*****************************************************************************************
 * @file: Rotation.cpp
 *
 * @brief: This file provides rotation matrices in the form of a templated Rotation class.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Rotation.hpp"

#include <cmath>
#include <stdexcept>

namespace sim
{

Rotation::Rotation()
{
    m_matrix = eye(3);
};

Rotation::Rotation(real_t roll, real_t pitch, real_t yaw)
{
    m_matrix = (Rz(yaw) * Ry(pitch) * Rx(roll)).m_matrix;
}

Rotation::Rotation(EulerAngles angles)
{
    m_matrix = (Rz(angles.yaw) * Ry(angles.pitch) * Rx(angles.roll)).m_matrix;
}

Rotation::Rotation(const Matrix& mat)
{
    dimension_t expected = {3, 3};
    if (size(mat) != expected)
    {
        throw std::invalid_argument("Incompatible dimension for rotation matrix");
    }
    m_matrix = mat;
}

Matrix Rotation::getMatrix() const
{
    return m_matrix;
}

EulerAngles Rotation::eulerAngles() const
{
    real_t sy = -m_matrix(3, 1);
    real_t cy = std::sqrt(m_matrix(1, 1) * m_matrix(1, 1) + m_matrix(2, 1) * m_matrix(2, 1));
    EulerAngles angles;
    angles.yaw = std::atan2(m_matrix(2, 1), m_matrix(1, 1));
    angles.pitch = std::atan2(sy, cy);
    angles.roll = std::atan2(m_matrix(3, 2), m_matrix(3, 3));
    return angles;
}

std::array<Vector, 3UL> Rotation::unitVectors() const
{
    std::array<Vector, 3UL> res;
    for (index_t c = 1; c <= static_cast<index_t>(3); ++c)
    {
        Vector col(3);
        for (index_t r = 1; r <= static_cast<index_t>(3); ++r)
        {
            col(r) = m_matrix(r, c);
        }
        res[static_cast<uint_t>(c - 1)] = (col);
    }
    return res;
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

Rotation Rotation::Rx(real_t roll)
{
    return Rotation(Matrix{{1, 0, 0}, {0, std::cos(roll), -std::sin(roll)}, {0, std::sin(roll), std::cos(roll)}});
}

Rotation Rotation::Ry(real_t pitch)
{
    return Rotation(Matrix{{std::cos(pitch), 0, std::sin(pitch)}, {0, 1, 0}, {-std::sin(pitch), 0, std::cos(pitch)}});
}

Rotation Rotation::Rz(real_t yaw)
{
    return Rotation(Matrix{{std::cos(yaw), -std::sin(yaw), 0}, {std::sin(yaw), std::cos(yaw), 0}, {0, 0, 1}});
}

std::ostream& operator<<(std::ostream& os, const Rotation& R)
{
    os << R.m_matrix;
    return os;
}

}  // namespace sim