/*****************************************************************************************
 * @file: Shapes.cpp
 *
 * @brief: This file implements the Shape class functions.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Shapes.hpp"

#include <cmath>
#include <stdexcept>

#include "Vector.hpp"

namespace sim
{

/*****************************************************************************************
 * Sphere Implementation
 ****************************************************************************************/

Sphere::Sphere(real_t radius) : m_radius(radius)
{
  if (radius <= 0) {
    throw std::invalid_argument("Radius must be positive");
  }
}

real_t Sphere::radius() const
{
  return m_radius;
}

real_t Sphere::volume() const
{
  return (4.0 / 3.0) * M_PI * m_radius * m_radius * m_radius;
}

Vector Sphere::centroid() const
{
  return Vector({0.0, 0.0, 0.0});
}

real_t Sphere::inertia(real_t rho) const
{
  if (rho <= 0) {
    throw std::invalid_argument("Density must be positive");
  }
  real_t m = rho * volume();
  return (2.0 / 5.0) * m * m_radius * m_radius;
}

std::vector<Vector> Sphere::vertices(real_t delta) const
{
  std::vector<Vector> verts;
  if (delta <= 0)
    delta = M_PI / 8.0;
  size_t n_phi = static_cast<size_t>(std::ceil(2.0 * M_PI / delta));
  size_t n_theta = static_cast<size_t>(std::ceil(M_PI / delta));
  real_t dphi = 2.0 * M_PI / static_cast<real_t>(n_phi);
  real_t dtheta = M_PI / static_cast<real_t>(n_theta);

  for (size_t ti = 0; ti <= n_theta; ++ti) {
    real_t theta = static_cast<real_t>(ti) * dtheta;  // 0..PI
    for (size_t pi = 0; pi < n_phi; ++pi) {
      real_t phi = static_cast<real_t>(pi) * dphi;  // 0..2PI
      real_t x = m_radius * std::sin(theta) * std::cos(phi);
      real_t y = m_radius * std::sin(theta) * std::sin(phi);
      real_t z = m_radius * std::cos(theta);
      verts.push_back(Vector({x, y, z}));
    }
  }
  return verts;
}

std::vector<Vector> Sphere::vertices() const
{
  return vertices(M_PI / 8.0);
}

bool_t Sphere::containsPoint(const Vector& point) const
{
  return norm(point) <= m_radius;
}

/*****************************************************************************************
 * Box (cuboid) Implementation
 ****************************************************************************************/

Box::Box(real_t width, real_t height, real_t depth) : m_width(width), m_height(height), m_depth(depth)
{
  if (width <= 0) {
    throw std::invalid_argument("Width must be positive");
  }
  if (height <= 0) {
    throw std::invalid_argument("Height must be positive");
  }
  if (depth <= 0) {
    throw std::invalid_argument("Depth must be positive");
  }
}

Box::Box(real_t size) : Box(size, size, size) {}

real_t Box::width() const
{
  return m_width;
}
real_t Box::height() const
{
  return m_height;
}
real_t Box::depth() const
{
  return m_depth;
}

real_t Box::volume() const
{
  return m_width * m_height * m_depth;
}

Vector Box::centroid() const
{
  return Vector({0.0, 0.0, 0.0});
}

real_t Box::inertia(real_t rho) const
{
  if (rho <= 0) {
    throw std::invalid_argument("Density must be positive");
  }
  real_t m = rho * volume();
  real_t Ix = (1.0 / 12.0) * m * (m_height * m_height + m_depth * m_depth);
  real_t Iy = (1.0 / 12.0) * m * (m_width * m_width + m_depth * m_depth);
  real_t Iz = (1.0 / 12.0) * m * (m_width * m_width + m_height * m_height);
  return (Ix + Iy + Iz) / 3.0;  // representative scalar
}

std::vector<Vector> Box::vertices() const
{
  real_t hw = m_width / 2.0;
  real_t hh = m_height / 2.0;
  real_t hd = m_depth / 2.0;
  return {
      Vector({-hw, -hh, -hd}),
      Vector({hw, -hh, -hd}),
      Vector({hw, hh, -hd}),
      Vector({-hw, hh, -hd}),  // bottom face
      Vector({-hw, -hh, hd}),
      Vector({hw, -hh, hd}),
      Vector({hw, hh, hd}),
      Vector({-hw, hh, hd})  // top face
  };
}

bool_t Box::containsPoint(const Vector& point) const
{
  return std::abs(point[0]) <= m_width / 2.0 && std::abs(point[1]) <= m_height / 2.0 &&
         std::abs(point[2]) <= m_depth / 2.0;
}

}  // namespace sim
