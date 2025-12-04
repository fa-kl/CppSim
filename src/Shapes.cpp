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
 * Circle Implementation
 ****************************************************************************************/

Circle::Circle(real_t radius) : m_radius(radius)
{
  if (radius <= 0) {
    throw std::invalid_argument("Radius must be positive");
  }
}

real_t Circle::radius() const
{
  return m_radius;
}

real_t Circle::area() const
{
  return M_PI * m_radius * m_radius;
}

Vector Circle::centroid() const
{
  return Vector({0.0, 0.0});
}

real_t Circle::inertia(real_t rho) const
{
  if (rho <= 0) {
    throw std::invalid_argument("Density must be positive");
  }
  real_t m = rho * area();
  return 0.5 * m * m_radius * m_radius;
}

std::vector<Vector> Circle::vertices(real_t delta_phi) const
{
  std::vector<Vector> verts;
  size_t n = static_cast<size_t>(std::ceil(2.0 * M_PI / delta_phi));
  delta_phi = 2.0 * M_PI / static_cast<real_t>(n);  // Adjust to get exact 2π coverage

  for (size_t i = 0; i < n; ++i) {
    real_t phi = static_cast<real_t>(i) * delta_phi;
    verts.push_back(Vector({m_radius * std::cos(phi), m_radius * std::sin(phi)}));
  }
  return verts;
}

std::vector<Vector> Circle::vertices() const
{
  return vertices(M_PI / 8.0);
}

bool_t Circle::containsPoint(const Vector& point) const
{
  return norm(point) <= m_radius;
}

/*****************************************************************************************
 * Box Implementation
 ****************************************************************************************/

Box::Box(real_t width, real_t height) : m_width(width), m_height(height)
{
  if (width <= 0) {
    throw std::invalid_argument("Width must be positive");
  }
  if (height <= 0) {
    throw std::invalid_argument("Height must be positive");
  }
}

Box::Box(real_t size) : Box(size, size) {}

real_t Box::width() const
{
  return m_width;
}

real_t Box::height() const
{
  return m_height;
}

real_t Box::area() const
{
  return m_width * m_height;
}

Vector Box::centroid() const
{
  return Vector({0.0, 0.0});
}

real_t Box::inertia(real_t rho) const
{
  if (rho <= 0) {
    throw std::invalid_argument("Density must be positive");
  }
  real_t m = rho * area();
  return (m / 12.0) * (m_width * m_width + m_height * m_height);
}

std::vector<Vector> Box::vertices() const
{
  real_t hw = m_width / 2.0;
  real_t hh = m_height / 2.0;
  return {Vector({-hw, -hh}),  // Bottom-left
          Vector({hw, -hh}),   // Bottom-right
          Vector({hw, hh}),    // Top-right
          Vector({-hw, hh})};  // Top-left
}

bool_t Box::containsPoint(const Vector& point) const
{
  return std::abs(point[0]) <= m_width / 2.0 && std::abs(point[1]) <= m_height / 2.0;
}

/*****************************************************************************************
 * Polygon Implementation
 ****************************************************************************************/

Polygon::Polygon(const std::vector<Vector>& vertices) : m_vertices(vertices)
{
  if (vertices.size() < 3) {
    throw std::invalid_argument("Polygon must have at least 3 vertices");
  }
  // Verify all vertices are 2D
  for (const auto& v : vertices) {
    if (v.length() != 2) {
      throw std::invalid_argument("All vertices must be 2-dimensional");
    }
  }
}

Polygon::Polygon(std::initializer_list<Vector> vertices) : Polygon(std::vector<Vector>(vertices)) {}

const std::vector<Vector>& Polygon::getVertices() const
{
  return m_vertices;
}

real_t Polygon::area() const
{
  size_t n = m_vertices.size();
  real_t A = 0.0;
  for (size_t i = 0; i < n; ++i) {
    size_t j = (i + 1) % n;
    A += m_vertices[i][0] * m_vertices[j][1];
    A -= m_vertices[j][0] * m_vertices[i][1];
  }
  return std::abs(A) / 2.0;
}

Vector Polygon::centroid() const
{
  size_t n = m_vertices.size();
  real_t A = area();
  real_t cx = 0.0, cy = 0.0;

  for (size_t i = 0; i < n; ++i) {
    size_t j = (i + 1) % n;
    real_t cross = m_vertices[i][0] * m_vertices[j][1] - m_vertices[j][0] * m_vertices[i][1];
    cx += (m_vertices[i][0] + m_vertices[j][0]) * cross;
    cy += (m_vertices[i][1] + m_vertices[j][1]) * cross;
  }

  return Vector({cx / (6.0 * A), cy / (6.0 * A)});
}

real_t Polygon::inertia(real_t rho) const
{
  if (rho <= 0) {
    throw std::invalid_argument("Density must be positive");
  }

  size_t n = m_vertices.size();
  Vector c = centroid();

  // Translate vertices to centroid frame
  std::vector<Vector> v_centered;
  v_centered.reserve(n);
  for (const auto& v : m_vertices) {
    v_centered.push_back(v - c);
  }

  // Compute inertia about centroid
  real_t I = 0.0;
  for (size_t i = 0; i < n; ++i) {
    size_t j = (i + 1) % n;
    const Vector& v1 = v_centered[i];
    const Vector& v2 = v_centered[j];

    real_t cross = v1[0] * v2[1] - v2[0] * v1[1];
    real_t term = v1[0] * v1[0] + v1[0] * v2[0] + v2[0] * v2[0] + v1[1] * v1[1] + v1[1] * v2[1] + v2[1] * v2[1];
    I += cross * term;
  }

  return rho * std::abs(I) / 12.0;
}

std::vector<Vector> Polygon::vertices() const
{
  return m_vertices;
}

bool_t Polygon::containsPoint(const Vector& point) const
{
  // Winding number algorithm
  int winding = 0;
  size_t n = m_vertices.size();

  for (size_t i = 0; i < n; ++i) {
    size_t j = (i + 1) % n;
    const Vector& v1 = m_vertices[i];
    const Vector& v2 = m_vertices[j];

    if (v1[1] <= point[1]) {
      if (v2[1] > point[1]) {  // Upward crossing
        real_t cross = (v2[0] - v1[0]) * (point[1] - v1[1]) - (point[0] - v1[0]) * (v2[1] - v1[1]);
        if (cross > 0) {
          ++winding;
        }
      }
    } else {
      if (v2[1] <= point[1]) {  // Downward crossing
        real_t cross = (v2[0] - v1[0]) * (point[1] - v1[1]) - (point[0] - v1[0]) * (v2[1] - v1[1]);
        if (cross < 0) {
          --winding;
        }
      }
    }
  }

  return winding != 0;
}

}  // namespace sim
