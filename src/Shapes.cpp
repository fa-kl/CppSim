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

Sphere::Sphere(real_t radius, Color color) : m_radius(radius), m_color(color)
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

std::vector<Vertex> Sphere::vertices(real_t delta) const
{
  std::vector<Vertex> verts;
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
      verts.push_back(Vertex(Vector({x, y, z}), m_color));
    }
  }
  return verts;
}

std::vector<Vertex> Sphere::vertices() const
{
  return vertices(M_PI / 8.0);
}

bool_t Sphere::containsPoint(const Vector& point) const
{
  return norm(point) <= m_radius;
}

Mesh Sphere::mesh() const
{
  Mesh result;
  std::vector<Vertex> verts = vertices();
  if (verts.size() < 3) {
    return result;
  }
  // Create triangles using fan triangulation from first vertex
  for (size_t i = 1; i + 1 < verts.size(); ++i) {
    result.push_back(Triangle{verts[0], verts[i], verts[i + 1]});
  }
  return result;
}

/*****************************************************************************************
 * Box (cuboid) Implementation
 ****************************************************************************************/

Box::Box(real_t width, real_t height, real_t depth, Color color)
    : m_width(width), m_height(height), m_depth(depth), m_color(color)
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

Box::Box(real_t size, Color color) : Box(size, size, size, color) {}

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

std::vector<Vertex> Box::vertices() const
{
  real_t hw = m_width / 2.0;
  real_t hh = m_height / 2.0;
  real_t hd = m_depth / 2.0;
  return {Vertex(Vector({-hw, -hh, -hd}), m_color),
          Vertex(Vector({hw, -hh, -hd}), m_color),
          Vertex(Vector({hw, hh, -hd}), m_color),
          Vertex(Vector({-hw, hh, -hd}), m_color),
          Vertex(Vector({-hw, -hh, hd}), m_color),
          Vertex(Vector({hw, -hh, hd}), m_color),
          Vertex(Vector({hw, hh, hd}), m_color),
          Vertex(Vector({-hw, hh, hd}), m_color)};
}

Mesh Box::mesh() const
{
  Mesh result;
  std::vector<Vertex> verts = vertices();
  result.push_back(Triangle{verts[0], verts[1], verts[2]});
  result.push_back(Triangle{verts[0], verts[2], verts[3]});
  result.push_back(Triangle{verts[4], verts[6], verts[5]});
  result.push_back(Triangle{verts[4], verts[7], verts[6]});
  result.push_back(Triangle{verts[0], verts[4], verts[5]});
  result.push_back(Triangle{verts[0], verts[5], verts[1]});
  result.push_back(Triangle{verts[2], verts[6], verts[7]});
  result.push_back(Triangle{verts[2], verts[7], verts[3]});
  result.push_back(Triangle{verts[0], verts[3], verts[7]});
  result.push_back(Triangle{verts[0], verts[7], verts[4]});
  result.push_back(Triangle{verts[1], verts[5], verts[6]});
  result.push_back(Triangle{verts[1], verts[6], verts[2]});
  return result;
}

bool_t Box::containsPoint(const Vector& point) const
{
  return std::abs(point[0]) <= m_width / 2.0 && std::abs(point[1]) <= m_height / 2.0 &&
         std::abs(point[2]) <= m_depth / 2.0;
}

}  // namespace sim
