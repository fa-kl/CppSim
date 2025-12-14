/*****************************************************************************************
 * @file: RigidBody.hpp
 *
 * @brief: This file implements a rigid body class.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "RigidBody.hpp"

#include <stdexcept>

namespace sim
{

RigidBody::RigidBody(Mesh mesh, Material material, Transform transform)
    : m_mesh(mesh),
      m_mass(mesh.mass(material.density)),
      m_inertia(mesh.inertia(material.density)),
      m_transform(transform),
      m_material(material)
{
    if (m_mass <= 0)
    {
        throw std::invalid_argument(
            "Well this is unexpected. The mass of a rigid body cannot be negative, but the internal computations "
            "resulted in a negative mass.");
    }
}

const Mesh& RigidBody::getMesh() const
{
    return m_mesh;
}

const Transform& RigidBody::getTransform() const
{
    return m_transform;
}

void RigidBody::setTransform(const Transform& transform)
{
    m_transform = transform;
}

const Material& RigidBody::getMaterial() const
{
    return m_material;
}

real_t RigidBody::getMass() const
{
    return m_mass;
}

const Matrix& RigidBody::getInertia() const
{
    return m_inertia;
}

}  // namespace sim