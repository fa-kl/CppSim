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
    if (radius <= 0)
    {
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
    if (rho <= 0)
    {
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

    for (size_t ti = 0; ti <= n_theta; ++ti)
    {
        real_t theta = static_cast<real_t>(ti) * dtheta;  // 0..PI
        for (size_t pi = 0; pi < n_phi; ++pi)
        {
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

Mesh Sphere::mesh(real_t delta) const
{
    Mesh result;
    std::vector<Vertex> verts = vertices(delta);
    if (verts.size() < 3)
    {
        return result;
    }
    if (delta <= 0)
    {
        delta = M_PI / 16.0;
    }
    size_t n_phi = static_cast<size_t>(std::ceil(2.0 * M_PI / delta));
    size_t n_theta = static_cast<size_t>(std::ceil(M_PI / delta));
    for (size_t ti = 0; ti < n_theta; ++ti)
    {
        for (size_t pi = 0; pi < n_phi; ++pi)
        {
            size_t i0 = ti * n_phi + pi;
            size_t i1 = ti * n_phi + ((pi + 1) % n_phi);
            size_t i2 = (ti + 1) * n_phi + pi;
            size_t i3 = (ti + 1) * n_phi + ((pi + 1) % n_phi);
            result.push_back(Triangle{verts[i0], verts[i1], verts[i2]});
            result.push_back(Triangle{verts[i1], verts[i3], verts[i2]});
        }
    }

    return result;
}

Mesh Sphere::mesh() const
{
    return mesh(M_PI / 16.0);
}

/*****************************************************************************************
 * Box (cuboid) Implementation
 ****************************************************************************************/

Box::Box(real_t width, real_t height, real_t depth, Color color)
    : m_width(width), m_height(height), m_depth(depth), m_color(color)
{
    if (width <= 0)
    {
        throw std::invalid_argument("Width must be positive");
    }
    if (height <= 0)
    {
        throw std::invalid_argument("Height must be positive");
    }
    if (depth <= 0)
    {
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
    if (rho <= 0)
    {
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

/*****************************************************************************************
 * Cylinder Implementation
 ****************************************************************************************/

Cylinder::Cylinder(real_t radius, real_t height, Color color) : m_radius(radius), m_height(height), m_color(color)
{
    if (radius <= 0)
    {
        throw std::invalid_argument("Radius must be positive");
    }
    if (height <= 0)
    {
        throw std::invalid_argument("Height must be positive");
    }
}

real_t Cylinder::radius() const
{
    return m_radius;
}

real_t Cylinder::height() const
{
    return m_height;
}

real_t Cylinder::volume() const
{
    return M_PI * m_radius * m_radius * m_height;
}

Vector Cylinder::centroid() const
{
    return Vector({0.0, 0.0, 0.0});
}

real_t Cylinder::inertia(real_t rho) const
{
    if (rho <= 0)
    {
        throw std::invalid_argument("Density must be positive");
    }
    real_t m = rho * volume();
    real_t Iz = 0.5 * m * m_radius * m_radius;
    real_t Ix = (1.0 / 12.0) * m * (3.0 * m_radius * m_radius + m_height * m_height);
    return (Iz + 2.0 * Ix) / 3.0;
}

std::vector<Vertex> Cylinder::vertices(real_t delta) const
{
    std::vector<Vertex> verts;
    if (delta <= 0)
        delta = M_PI / 8.0;
    size_t n = static_cast<size_t>(std::ceil(2.0 * M_PI / delta));
    real_t dtheta = 2.0 * M_PI / static_cast<real_t>(n);
    real_t half_h = m_height / 2.0;

    // Bottom circle
    for (size_t i = 0; i < n; ++i)
    {
        real_t theta = static_cast<real_t>(i) * dtheta;
        real_t x = m_radius * std::cos(theta);
        real_t y = m_radius * std::sin(theta);
        verts.push_back(Vertex(Vector({x, y, -half_h}), m_color));
    }

    // Top circle
    for (size_t i = 0; i < n; ++i)
    {
        real_t theta = static_cast<real_t>(i) * dtheta;
        real_t x = m_radius * std::cos(theta);
        real_t y = m_radius * std::sin(theta);
        verts.push_back(Vertex(Vector({x, y, half_h}), m_color));
    }

    // Center points for caps
    verts.push_back(Vertex(Vector({0.0, 0.0, -half_h}), m_color));  // bottom center
    verts.push_back(Vertex(Vector({0.0, 0.0, half_h}), m_color));   // top center

    return verts;
}

std::vector<Vertex> Cylinder::vertices() const
{
    return vertices(M_PI / 8.0);
}

Mesh Cylinder::mesh(real_t delta) const
{
    Mesh result;
    std::vector<Vertex> verts = vertices(delta);
    if (delta <= 0)
        delta = M_PI / 8.0;
    size_t n = static_cast<size_t>(std::ceil(2.0 * M_PI / delta));

    size_t bottom_center = 2 * n;
    size_t top_center = 2 * n + 1;

    // Side triangles
    for (size_t i = 0; i < n; ++i)
    {
        size_t next = (i + 1) % n;
        result.push_back(Triangle{verts[i], verts[i + n], verts[next]});
        result.push_back(Triangle{verts[next], verts[i + n], verts[next + n]});
    }

    // Bottom cap
    for (size_t i = 0; i < n; ++i)
    {
        size_t next = (i + 1) % n;
        result.push_back(Triangle{verts[bottom_center], verts[next], verts[i]});
    }

    // Top cap
    for (size_t i = 0; i < n; ++i)
    {
        size_t next = (i + 1) % n;
        result.push_back(Triangle{verts[top_center], verts[i + n], verts[next + n]});
    }

    return result;
}

Mesh Cylinder::mesh() const
{
    return mesh(M_PI / 8.0);
}

bool_t Cylinder::containsPoint(const Vector& point) const
{
    real_t r = std::sqrt(point[0] * point[0] + point[1] * point[1]);
    return r <= m_radius && std::abs(point[2]) <= m_height / 2.0;
}

/*****************************************************************************************
 * Pyramid Implementation
 ****************************************************************************************/

Pyramid::Pyramid(real_t base, real_t height, Color color) : m_base(base), m_height(height), m_color(color)
{
    if (base <= 0)
    {
        throw std::invalid_argument("Base must be positive");
    }
    if (height <= 0)
    {
        throw std::invalid_argument("Height must be positive");
    }
}

real_t Pyramid::base() const
{
    return m_base;
}

real_t Pyramid::height() const
{
    return m_height;
}

real_t Pyramid::volume() const
{
    return (1.0 / 3.0) * m_base * m_base * m_height;
}

Vector Pyramid::centroid() const
{
    return Vector({0.0, 0.0, -m_height / 4.0});  // Centroid is 1/4 up from base
}

real_t Pyramid::inertia(real_t rho) const
{
    if (rho <= 0)
    {
        throw std::invalid_argument("Density must be positive");
    }
    real_t m = rho * volume();
    // Approximate with simplified formula
    return m * (m_base * m_base + m_height * m_height) / 20.0;
}

std::vector<Vertex> Pyramid::vertices() const
{
    real_t hb = m_base / 2.0;
    real_t half_h = m_height / 2.0;
    return {Vertex(Vector({-hb, -hb, -half_h}), m_color),  // Base corners
            Vertex(Vector({hb, -hb, -half_h}), m_color),
            Vertex(Vector({hb, hb, -half_h}), m_color),
            Vertex(Vector({-hb, hb, -half_h}), m_color),
            Vertex(Vector({0.0, 0.0, half_h}), m_color)};  // Apex
}

Mesh Pyramid::mesh() const
{
    Mesh result;
    std::vector<Vertex> verts = vertices();

    // Base
    result.push_back(Triangle{verts[0], verts[1], verts[2]});
    result.push_back(Triangle{verts[0], verts[2], verts[3]});

    // Sides
    result.push_back(Triangle{verts[0], verts[4], verts[1]});
    result.push_back(Triangle{verts[1], verts[4], verts[2]});
    result.push_back(Triangle{verts[2], verts[4], verts[3]});
    result.push_back(Triangle{verts[3], verts[4], verts[0]});

    return result;
}

bool_t Pyramid::containsPoint(const Vector& point) const
{
    real_t half_h = m_height / 2.0;
    if (point[2] < -half_h || point[2] > half_h)
        return false;

    // Scale factor at this height
    real_t t = (half_h - point[2]) / m_height;
    real_t max_xy = t * m_base / 2.0;

    return std::abs(point[0]) <= max_xy && std::abs(point[1]) <= max_xy;
}

/*****************************************************************************************
 * Cone Implementation
 ****************************************************************************************/

Cone::Cone(real_t radius, real_t height, Color color) : m_radius(radius), m_height(height), m_color(color)
{
    if (radius <= 0)
    {
        throw std::invalid_argument("Radius must be positive");
    }
    if (height <= 0)
    {
        throw std::invalid_argument("Height must be positive");
    }
}

real_t Cone::radius() const
{
    return m_radius;
}

real_t Cone::height() const
{
    return m_height;
}

real_t Cone::volume() const
{
    return (1.0 / 3.0) * M_PI * m_radius * m_radius * m_height;
}

Vector Cone::centroid() const
{
    return Vector({0.0, 0.0, -m_height / 4.0});
}

real_t Cone::inertia(real_t rho) const
{
    if (rho <= 0)
    {
        throw std::invalid_argument("Density must be positive");
    }
    real_t m = rho * volume();
    real_t Iz = (3.0 / 10.0) * m * m_radius * m_radius;
    real_t Ix = (3.0 / 20.0) * m * (m_radius * m_radius + 2.0 * m_height * m_height);
    return (Iz + 2.0 * Ix) / 3.0;
}

std::vector<Vertex> Cone::vertices(real_t delta) const
{
    std::vector<Vertex> verts;
    if (delta <= 0)
        delta = M_PI / 8.0;
    size_t n = static_cast<size_t>(std::ceil(2.0 * M_PI / delta));
    real_t dtheta = 2.0 * M_PI / static_cast<real_t>(n);
    real_t half_h = m_height / 2.0;

    // Base circle
    for (size_t i = 0; i < n; ++i)
    {
        real_t theta = static_cast<real_t>(i) * dtheta;
        real_t x = m_radius * std::cos(theta);
        real_t y = m_radius * std::sin(theta);
        verts.push_back(Vertex(Vector({x, y, -half_h}), m_color));
    }

    // Apex
    verts.push_back(Vertex(Vector({0.0, 0.0, half_h}), m_color));

    // Base center
    verts.push_back(Vertex(Vector({0.0, 0.0, -half_h}), m_color));

    return verts;
}

std::vector<Vertex> Cone::vertices() const
{
    return vertices(M_PI / 8.0);
}

Mesh Cone::mesh(real_t delta) const
{
    Mesh result;
    std::vector<Vertex> verts = vertices(delta);
    if (delta <= 0)
        delta = M_PI / 8.0;
    size_t n = static_cast<size_t>(std::ceil(2.0 * M_PI / delta));

    size_t apex = n;
    size_t base_center = n + 1;

    // Side triangles
    for (size_t i = 0; i < n; ++i)
    {
        size_t next = (i + 1) % n;
        result.push_back(Triangle{verts[i], verts[apex], verts[next]});
    }

    // Base cap
    for (size_t i = 0; i < n; ++i)
    {
        size_t next = (i + 1) % n;
        result.push_back(Triangle{verts[base_center], verts[next], verts[i]});
    }

    return result;
}

Mesh Cone::mesh() const
{
    return mesh(M_PI / 8.0);
}

bool_t Cone::containsPoint(const Vector& point) const
{
    real_t half_h = m_height / 2.0;
    if (point[2] < -half_h || point[2] > half_h)
        return false;

    real_t t = (half_h - point[2]) / m_height;
    real_t max_r = t * m_radius;
    real_t r = std::sqrt(point[0] * point[0] + point[1] * point[1]);

    return r <= max_r;
}

/*****************************************************************************************
 * Polytope Implementation
 ****************************************************************************************/

Polytope::Polytope(const std::vector<Vector>& vertices, Color color) : m_color(color), m_centroid_computed(false)
{
    if (vertices.size() < 4)
    {
        throw std::invalid_argument("Polytope requires at least 4 vertices");
    }
    for (const auto& v : vertices)
    {
        if (v.length() != 3)
        {
            throw std::invalid_argument("All vertices must be 3-dimensional");
        }
        m_vertices.push_back(v);
    }
}

const std::vector<Vector>& Polytope::getVertices() const
{
    return m_vertices;
}

real_t Polytope::volume() const
{
    // Simplified volume calculation: divide into tetrahedra from centroid
    Vector c = centroid();
    real_t vol = 0.0;

    // This is a simplified approximation - proper convex hull needed for accuracy
    Mesh m = mesh();
    for (const auto& tri : m)
    {
        Vector v0 = tri[0].point - c;
        Vector v1 = tri[1].point - c;
        Vector v2 = tri[2].point - c;
        vol += std::abs(dot(v0, cross(v1, v2))) / 6.0;
    }

    return vol;
}

Vector Polytope::centroid() const
{
    if (!m_centroid_computed)
    {
        m_centroid = Vector({0.0, 0.0, 0.0});
        for (const auto& v : m_vertices)
        {
            m_centroid = m_centroid + v;
        }
        m_centroid = m_centroid / static_cast<real_t>(m_vertices.size());
        m_centroid_computed = true;
    }
    return m_centroid;
}

real_t Polytope::inertia(real_t rho) const
{
    if (rho <= 0)
    {
        throw std::invalid_argument("Density must be positive");
    }
    // Simplified approximation
    real_t m = rho * volume();
    Vector c = centroid();
    real_t I = 0.0;
    for (const auto& v : m_vertices)
    {
        Vector r = v - c;
        I += norm(r) * norm(r);
    }
    return m * I / static_cast<real_t>(m_vertices.size());
}

std::vector<Vertex> Polytope::vertices() const
{
    std::vector<Vertex> verts;
    for (const auto& v : m_vertices)
    {
        verts.push_back(Vertex(v, m_color));
    }
    return verts;
}

Mesh Polytope::mesh() const
{
    // Simple convex hull approximation using fan triangulation from centroid
    // For a proper implementation, use a convex hull algorithm (e.g., QuickHull)
    Mesh result;
    std::vector<Vertex> verts = vertices();

    if (verts.size() < 4)
        return result;

    // Simplified: create triangles using first 3 vertices as base and fan from others
    for (size_t i = 3; i < verts.size(); ++i)
    {
        result.push_back(Triangle{verts[0], verts[1], verts[i]});
        result.push_back(Triangle{verts[0], verts[i], verts[2]});
        result.push_back(Triangle{verts[1], verts[2], verts[i]});
    }

    return result;
}

bool_t Polytope::containsPoint(const Vector& point) const
{
    // Simplified point-in-polytope test
    // Proper implementation would test against all faces
    Vector c = centroid();
    real_t max_dist = 0.0;
    for (const auto& v : m_vertices)
    {
        max_dist = std::max(max_dist, norm(v - c));
    }
    return norm(point - c) <= max_dist;
}

}  // namespace sim
