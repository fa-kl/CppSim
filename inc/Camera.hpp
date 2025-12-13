/*****************************************************************************************
 * @file: Camera.hpp
 *
 * @brief: Camera class for 3D rendering and projection.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include "Matrix.hpp"
#include "Transform.hpp"
#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

class Camera
{
  protected:
    real_t m_fov;
    real_t m_near;
    real_t m_far;
    real_t m_aspect;
    Transform m_transform;
    Matrix m_projection;

  public:
    Camera(const Transform& transform,
           real_t aspect = 16.0 / 9.0,
           real_t fov = deg2rad(60),
           real_t near = 0.1,
           real_t far = 100);

    Vector projectWorldToCamera(const Vector& world_point) const;

    Vector projectWorldToNDC(const Vector& world_point) const;

    bool isInView(const Vector& world_point) const;

    Transform getTransform() const;

    void setTransform(const Transform& transform);

    Matrix getProjectionMatrix() const;

    real_t getFOV() const;

    real_t getAspect() const;

    real_t getNear() const;

    real_t getFar() const;
};

}  // namespace sim
