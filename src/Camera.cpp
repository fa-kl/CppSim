/*****************************************************************************************
 * @file: Camera.cpp
 *
 * @brief: Camera implementation.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Camera.hpp"

#include <cmath>
#include <stdexcept>

namespace sim
{

Camera::Camera(const Transform& transform, real_t aspect, real_t fov, real_t near, real_t far)
    : m_fov(fov), m_near(near), m_far(far), m_aspect(aspect), m_transform(transform)
{
  if (aspect <= 0) {
    throw std::invalid_argument("Aspect ratio must be greater than zero");
  }
  if (fov <= 0) {
    throw std::invalid_argument("Field of view must be greater than zero");
  }
  if (near <= 0) {
    throw std::invalid_argument("Near distance must be greater than zero");
  }
  if (far <= 0) {
    throw std::invalid_argument("Far distance must be greater than zero");
  }
  // OpenGL-style perspective projection matrix
  real_t f = 1.0 / std::tan(m_fov / 2.0);
  real_t range = m_far - m_near;
  m_projection = Matrix(4, 4);
  m_projection(1, 1) = f / m_aspect;
  m_projection(2, 2) = f;
  m_projection(3, 3) = -(m_far + m_near) / range;
  m_projection(3, 4) = -2.0 * m_far * m_near / range;
  m_projection(4, 3) = -1.0;
  m_projection(4, 4) = 0.0;
}

Vector Camera::projectWorldToCamera(const Vector& world_point) const
{
  if (world_point.length() != 3) {
    throw std::invalid_argument("World point must be 3-dimensional");
  }
  Transform world_to_camera = inv(m_transform);
  Vector camera_pos = world_to_camera * world_point;
  return camera_pos;
}

Vector Camera::projectWorldToNDC(const Vector& world_point) const
{
  Vector camera_pos = projectWorldToCamera(world_point);
  Vector clip_space = m_projection * Vector({camera_pos(1), camera_pos(2), camera_pos(3), 1.0});
  if (std::abs(clip_space(4)) < 1e-10) {
    throw std::runtime_error("Division by zero in perspective divide");
  }
  Vector ndc = {clip_space(1) / clip_space(4), clip_space(2) / clip_space(4), clip_space(3) / clip_space(4)};
  return ndc;
}

bool Camera::isInView(const Vector& world_point) const
{
  try {
    Vector camera_pos = projectWorldToCamera(world_point);
    // Check if point is behind the camera (z > 0 in camera space, camera looks in -z direction)
    if (camera_pos(3) >= 0) {
      return false;
    }
    // Check if point is within near and far planes
    real_t depth = -camera_pos(3);  // Distance from camera
    if (depth < m_near || depth > m_far) {
      return false;
    }
    // Check if point is within the view frustum
    Vector ndc = projectWorldToNDC(world_point);
    // Point is in view if NDC coordinates are in [-1, 1] range
    if (ndc(1) < -1.0 || ndc(1) > 1.0 || ndc(2) < -1.0 || ndc(2) > 1.0 || ndc(3) < -1.0 || ndc(3) > 1.0) {
      return false;
    }
    return true;
  } catch (...) {
    return false;
  }
}

Transform Camera::getTransform() const
{
  return m_transform;
}

void Camera::setTransform(const Transform& transform)
{
  m_transform = transform;
}

Matrix Camera::getProjectionMatrix() const
{
  return m_projection;
}

real_t Camera::getFOV() const
{
  return m_fov;
}

real_t Camera::getAspect() const
{
  return m_aspect;
}

real_t Camera::getNear() const
{
  return m_near;
}

real_t Camera::getFar() const
{
  return m_far;
}

}  // namespace sim
