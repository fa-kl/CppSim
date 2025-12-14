/*****************************************************************************************
 * @file: RigidBody.hpp
 *
 * @brief: This file provides a rigid body class.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/
#pragma once

#include "Material.hpp"
#include "Matrix.hpp"
#include "Mesh.hpp"
#include "Transform.hpp"
#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

class RigidBody
{
  protected:
    Mesh m_mesh;
    real_t m_mass;
    Matrix m_inertia;
    Transform m_transform;
    Material m_material;

  public:
    RigidBody(Mesh mesh, Material material, Transform transform = Transform());

    /**
     * @brief Get the mesh of the rigid body.
     */
    const Mesh& getMesh() const;

    /**
     * @brief Get the transform of the rigid body.
     */
    const Transform& getTransform() const;

    /**
     * @brief Set the transform of the rigid body.
     */
    void setTransform(const Transform& transform);

    /**
     * @brief Get the material of the rigid body.
     */
    const Material& getMaterial() const;

    /**
     * @brief Get the mass of the rigid body.
     */
    real_t getMass() const;

    /**
     * @brief Get the inertia tensor of the rigid body.
     */
    const Matrix& getInertia() const;
};

}  // namespace sim