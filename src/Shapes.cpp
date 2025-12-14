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

/*================================================================================================
 * Box (cuboid) Implementation
 ================================================================================================*/

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

Vector Box::center_of_mass() const
{
    return centroid();
}

Matrix Box::inertia(real_t rho) const
{
    if (rho <= 0)
    {
        throw std::invalid_argument("Density must be positive");
    }
    const real_t m = rho * volume();
    real_t Ix = (1.0 / 12.0) * m * (m_height * m_height + m_depth * m_depth);
    real_t Iy = (1.0 / 12.0) * m * (m_width * m_width + m_depth * m_depth);
    real_t Iz = (1.0 / 12.0) * m * (m_width * m_width + m_height * m_height);
    return {{Ix, 0, 0}, {0, Iy, 0}, {0, 0, Iz}};
}

bool_t Box::containsPoint(const Vector& point) const
{
    return std::abs(point[0]) <= m_width / 2.0 && std::abs(point[1]) <= m_height / 2.0 &&
           std::abs(point[2]) <= m_depth / 2.0;
}

std::vector<Vertex> Box::vertices() const
{
    real_t hw = m_width / 2.0;
    real_t hh = m_height / 2.0;
    real_t hd = m_depth / 2.0;
    Vector neg_z = {0, 0, -1};  // Front face normal
    Vector pos_z = {0, 0, 1};   // Back face normal
    Vector neg_y = {0, -1, 0};  // Bottom face normal
    Vector pos_y = {0, 1, 0};   // Top face normal
    Vector neg_x = {-1, 0, 0};  // Left face normal
    Vector pos_x = {1, 0, 0};   // Right face normal

    return {Vertex(Vector({-hw, -hh, -hd}), m_color, neg_z),  // 0: Front-bottom-left
            Vertex(Vector({hw, -hh, -hd}), m_color, neg_z),   // 1: Front-bottom-right
            Vertex(Vector({hw, hh, -hd}), m_color, neg_z),    // 2: Front-top-right
            Vertex(Vector({-hw, hh, -hd}), m_color, neg_z),   // 3: Front-top-left
            Vertex(Vector({-hw, -hh, hd}), m_color, pos_z),   // 4: Back-bottom-left
            Vertex(Vector({hw, -hh, hd}), m_color, pos_z),    // 5: Back-bottom-right
            Vertex(Vector({hw, hh, hd}), m_color, pos_z),     // 6: Back-top-right
            Vertex(Vector({-hw, hh, hd}), m_color, pos_z)};   // 7: Back-top-left
}

std::vector<TriangleVertexIndices> Box::faces() const
{
    return {// Front face
            {0, 1, 2},
            {0, 2, 3},
            // Back face
            {4, 6, 5},
            {4, 7, 6},
            // Bottom face
            {0, 4, 5},
            {0, 5, 1},
            // Top face
            {2, 6, 7},
            {2, 7, 3},
            // Left face
            {0, 3, 7},
            {0, 7, 4},
            // Right face
            {1, 5, 6},
            {1, 6, 2}};
}

Mesh Box::mesh() const
{
    return Mesh(vertices(), faces());
}

/*================================================================================================
 * Sphere Implementation
 ================================================================================================*/

Sphere::Sphere(real_t radius, Color color, real_t delta) : m_radius(radius), m_color(color), m_delta(delta)
{
    if (radius <= 0)
    {
        throw std::invalid_argument("Radius must be positive");
    }
    if (delta <= 0)
    {
        m_delta = M_PI / 8.0;
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

Vector Sphere::center_of_mass() const
{
    return centroid();
}

Matrix Sphere::inertia(real_t rho) const
{
    if (rho <= 0)
    {
        throw std::invalid_argument("Density must be positive");
    }
    // For a sphere, the inertia tensor is isotropic: I = (2/5) * m * r^2
    real_t I_val = (2.0 / 5.0) * rho * volume() * m_radius * m_radius;

    // Create diagonal 3x3 inertia tensor
    Matrix I(3, 3);
    I(0, 0) = I_val;
    I(0, 1) = 0;
    I(0, 2) = 0;
    I(1, 0) = 0;
    I(1, 1) = I_val;
    I(1, 2) = 0;
    I(2, 0) = 0;
    I(2, 1) = 0;
    I(2, 2) = I_val;
    return I;
}

bool_t Sphere::containsPoint(const Vector& point) const
{
    return norm(point) <= m_radius;
}

std::vector<Vertex> Sphere::vertices() const
{
    std::vector<Vertex> verts;
    real_t delta = m_delta;
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
            Vector pos({x, y, z});
            // For a sphere, the normal at any point is the normalized position vector
            Vector normal = pos / m_radius;
            verts.push_back(Vertex(pos, m_color, normal));
        }
    }
    return verts;
}

std::vector<TriangleVertexIndices> Sphere::faces() const
{
    real_t delta = m_delta;
    size_t n_phi = static_cast<size_t>(std::ceil(2.0 * M_PI / delta));
    size_t n_theta = static_cast<size_t>(std::ceil(M_PI / delta));

    std::vector<TriangleVertexIndices> face_indices;
    face_indices.reserve(n_theta * n_phi * 2);

    for (size_t ti = 0; ti < n_theta; ++ti)
    {
        for (size_t pi = 0; pi < n_phi; ++pi)
        {
            size_t i0 = ti * n_phi + pi;
            size_t i1 = ti * n_phi + ((pi + 1) % n_phi);
            size_t i2 = (ti + 1) * n_phi + pi;
            size_t i3 = (ti + 1) * n_phi + ((pi + 1) % n_phi);
            face_indices.push_back({static_cast<uint_t>(i0), static_cast<uint_t>(i1), static_cast<uint_t>(i2)});
            face_indices.push_back({static_cast<uint_t>(i1), static_cast<uint_t>(i3), static_cast<uint_t>(i2)});
        }
    }

    return face_indices;
}

Mesh Sphere::mesh() const
{
    return Mesh(vertices(), faces());
}

/*================================================================================================
 * Cylinder Implementation
 ================================================================================================*/

Cylinder::Cylinder(real_t radius, real_t height, Color color, real_t delta)
    : m_radius(radius), m_height(height), m_color(color), m_delta(delta)
{
    if (radius <= 0)
    {
        throw std::invalid_argument("Radius must be positive");
    }
    if (height <= 0)
    {
        throw std::invalid_argument("Height must be positive");
    }
    if (delta <= 0)
    {
        m_delta = M_PI / 8.0;
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

Vector Cylinder::center_of_mass() const
{
    return centroid();
}

Matrix Cylinder::inertia(real_t rho) const
{
    if (rho <= 0)
    {
        throw std::invalid_argument("Density must be positive");
    }
    real_t m = rho * volume();
    real_t Iz = 0.5 * m * m_radius * m_radius;
    real_t Ix = (1.0 / 12.0) * m * (3.0 * m_radius * m_radius + m_height * m_height);
    real_t Iy = Ix;  // Same as Ix due to cylindrical symmetry

    // Create diagonal 3x3 inertia tensor
    Matrix I(3, 3);
    I(0, 0) = Ix;
    I(0, 1) = 0;
    I(0, 2) = 0;
    I(1, 0) = 0;
    I(1, 1) = Iy;
    I(1, 2) = 0;
    I(2, 0) = 0;
    I(2, 1) = 0;
    I(2, 2) = Iz;
    return I;
}

bool_t Cylinder::containsPoint(const Vector& point) const
{
    real_t r = std::sqrt(point[0] * point[0] + point[1] * point[1]);
    return r <= m_radius && std::abs(point[2]) <= m_height / 2.0;
}

std::vector<Vertex> Cylinder::vertices() const
{
    std::vector<Vertex> verts;
    real_t delta = m_delta;
    size_t n = static_cast<size_t>(std::ceil(2.0 * M_PI / delta));
    real_t dtheta = 2.0 * M_PI / static_cast<real_t>(n);
    real_t half_h = m_height / 2.0;

    // Bottom circle (radial normals pointing outward)
    for (size_t i = 0; i < n; ++i)
    {
        real_t theta = static_cast<real_t>(i) * dtheta;
        real_t x = m_radius * std::cos(theta);
        real_t y = m_radius * std::sin(theta);
        Vector pos({x, y, -half_h});
        Vector normal({std::cos(theta), std::sin(theta), 0});  // Radial normal
        verts.push_back(Vertex(pos, m_color, normal));
    }

    // Top circle (radial normals pointing outward)
    for (size_t i = 0; i < n; ++i)
    {
        real_t theta = static_cast<real_t>(i) * dtheta;
        real_t x = m_radius * std::cos(theta);
        real_t y = m_radius * std::sin(theta);
        Vector pos({x, y, half_h});
        Vector normal({std::cos(theta), std::sin(theta), 0});  // Radial normal
        verts.push_back(Vertex(pos, m_color, normal));
    }

    // Center points for caps
    Vector bottom_center_pos({0.0, 0.0, -half_h});
    verts.push_back(Vertex(bottom_center_pos, m_color, Vector({0, 0, -1})));  // bottom center, normal pointing down

    Vector top_center_pos({0.0, 0.0, half_h});
    verts.push_back(Vertex(top_center_pos, m_color, Vector({0, 0, 1})));  // top center, normal pointing up

    return verts;
}

std::vector<TriangleVertexIndices> Cylinder::faces() const
{
    real_t delta = m_delta;
    size_t n = static_cast<size_t>(std::ceil(2.0 * M_PI / delta));

    size_t bottom_center = 2 * n;
    size_t top_center = 2 * n + 1;

    std::vector<TriangleVertexIndices> face_indices;
    face_indices.reserve(n * 4);

    // Side triangles
    for (size_t i = 0; i < n; ++i)
    {
        size_t next = (i + 1) % n;
        face_indices.push_back({static_cast<uint_t>(i), static_cast<uint_t>(i + n), static_cast<uint_t>(next)});
        face_indices.push_back({static_cast<uint_t>(next), static_cast<uint_t>(i + n), static_cast<uint_t>(next + n)});
    }

    // Bottom cap
    for (size_t i = 0; i < n; ++i)
    {
        size_t next = (i + 1) % n;
        face_indices.push_back({static_cast<uint_t>(bottom_center), static_cast<uint_t>(next), static_cast<uint_t>(i)});
    }

    // Top cap
    for (size_t i = 0; i < n; ++i)
    {
        size_t next = (i + 1) % n;
        face_indices.push_back(
            {static_cast<uint_t>(top_center), static_cast<uint_t>(i + n), static_cast<uint_t>(next + n)});
    }

    return face_indices;
}

Mesh Cylinder::mesh() const
{
    return Mesh(vertices(), faces());
}

/*================================================================================================
 * Pyramid Implementation
 ================================================================================================*/

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

Vector Pyramid::center_of_mass() const
{
    return centroid();
}

Matrix Pyramid::inertia(real_t rho) const
{
    if (rho <= 0)
    {
        throw std::invalid_argument("Density must be positive");
    }
    real_t m = rho * volume();
    // Inertia tensor for a pyramid with base in xy-plane
    // Ixx = Iyy (symmetry) = (1/20) * m * (b^2 + h^2)
    // Izz = (1/10) * m * b^2
    real_t Ixy = (1.0 / 20.0) * m * (m_base * m_base + m_height * m_height);
    real_t Iz = (1.0 / 10.0) * m * m_base * m_base;

    // Create diagonal 3x3 inertia tensor
    Matrix I(3, 3);
    I(0, 0) = Ixy;
    I(0, 1) = 0;
    I(0, 2) = 0;
    I(1, 0) = 0;
    I(1, 1) = Ixy;
    I(1, 2) = 0;
    I(2, 0) = 0;
    I(2, 1) = 0;
    I(2, 2) = Iz;
    return I;
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

std::vector<Vertex> Pyramid::vertices() const
{
    real_t hb = m_base / 2.0;
    real_t half_h = m_height / 2.0;
    return {Vertex(Vector({-hb, -hb, -half_h}), m_color, Vector({0, 0, -1})),  // 0: Base corner - bottom normal
            Vertex(Vector({hb, -hb, -half_h}), m_color, Vector({0, 0, -1})),   // 1: Base corner - bottom normal
            Vertex(Vector({hb, hb, -half_h}), m_color, Vector({0, 0, -1})),    // 2: Base corner - bottom normal
            Vertex(Vector({-hb, hb, -half_h}), m_color, Vector({0, 0, -1})),   // 3: Base corner - bottom normal
            Vertex(Vector({0.0, 0.0, half_h}), m_color, Vector({0, 0, 1}))};   // 4: Apex - upward normal
}

std::vector<TriangleVertexIndices> Pyramid::faces() const
{
    return {// Base
            {0, 1, 2},
            {0, 2, 3},
            // Sides
            {0, 4, 1},
            {1, 4, 2},
            {2, 4, 3},
            {3, 4, 0}};
}

Mesh Pyramid::mesh() const
{
    return Mesh(vertices(), faces());
}

/*================================================================================================
 * Cone Implementation
 ================================================================================================*/

Cone::Cone(real_t radius, real_t height, Color color, real_t delta)
    : m_radius(radius), m_height(height), m_color(color), m_delta(delta)
{
    if (radius <= 0)
    {
        throw std::invalid_argument("Radius must be positive");
    }
    if (height <= 0)
    {
        throw std::invalid_argument("Height must be positive");
    }
    if (delta <= 0)
    {
        m_delta = M_PI / 8.0;
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

Vector Cone::center_of_mass() const
{
    return centroid();
}

Matrix Cone::inertia(real_t rho) const
{
    if (rho <= 0)
    {
        throw std::invalid_argument("Density must be positive");
    }
    real_t m = rho * volume();
    real_t Iz = (3.0 / 10.0) * m * m_radius * m_radius;
    real_t Ix = (3.0 / 20.0) * m * (m_radius * m_radius + 2.0 * m_height * m_height);
    real_t Iy = Ix;  // Same as Ix due to rotational symmetry

    // Create diagonal 3x3 inertia tensor
    Matrix I(3, 3);
    I(0, 0) = Ix;
    I(0, 1) = 0;
    I(0, 2) = 0;
    I(1, 0) = 0;
    I(1, 1) = Iy;
    I(1, 2) = 0;
    I(2, 0) = 0;
    I(2, 1) = 0;
    I(2, 2) = Iz;
    return I;
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

std::vector<Vertex> Cone::vertices() const
{
    std::vector<Vertex> verts;
    real_t delta = m_delta;
    size_t n = static_cast<size_t>(std::ceil(2.0 * M_PI / delta));
    real_t dtheta = 2.0 * M_PI / static_cast<real_t>(n);
    real_t half_h = m_height / 2.0;

    // Base circle (radial normals pointing outward)
    for (size_t i = 0; i < n; ++i)
    {
        real_t theta = static_cast<real_t>(i) * dtheta;
        real_t x = m_radius * std::cos(theta);
        real_t y = m_radius * std::sin(theta);
        Vector pos({x, y, -half_h});
        Vector normal({std::cos(theta), std::sin(theta), 0});  // Radial normal
        verts.push_back(Vertex(pos, m_color, normal));
    }

    // Apex
    verts.push_back(Vertex(Vector({0.0, 0.0, half_h}), m_color, Vector({0, 0, 1})));

    // Base center
    verts.push_back(Vertex(Vector({0.0, 0.0, -half_h}), m_color, Vector({0, 0, -1})));

    return verts;
}

std::vector<TriangleVertexIndices> Cone::faces() const
{
    real_t delta = m_delta;
    size_t n = static_cast<size_t>(std::ceil(2.0 * M_PI / delta));

    size_t apex = n;
    size_t base_center = n + 1;

    std::vector<TriangleVertexIndices> face_indices;
    face_indices.reserve(n * 2);

    // Side triangles
    for (size_t i = 0; i < n; ++i)
    {
        size_t next = (i + 1) % n;
        face_indices.push_back({static_cast<uint_t>(i), static_cast<uint_t>(apex), static_cast<uint_t>(next)});
    }

    // Base cap
    for (size_t i = 0; i < n; ++i)
    {
        size_t next = (i + 1) % n;
        face_indices.push_back({static_cast<uint_t>(base_center), static_cast<uint_t>(next), static_cast<uint_t>(i)});
    }

    return face_indices;
}

Mesh Cone::mesh() const
{
    return Mesh(vertices(), faces());
}

/*================================================================================================
 * Polytope Implementation
 ================================================================================================*/

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
    const auto& verts = m.getVertices();
    const auto& face_list = m.getFaceIndices();

    for (const auto& face : face_list)
    {
        Vector v0 = verts[face[0]].position - c;
        Vector v1 = verts[face[1]].position - c;
        Vector v2 = verts[face[2]].position - c;
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

Vector Polytope::center_of_mass() const
{
    return centroid();
}

Matrix Polytope::inertia(real_t rho) const
{
    if (rho <= 0)
    {
        throw std::invalid_argument("Density must be positive");
    }
    // Compute inertia tensor by averaging vertex contributions
    real_t m = rho * volume();
    Vector c = centroid();

    // Initialize 3x3 matrix to zero
    Matrix I(3, 3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            I(i, j) = 0.0;

    // Sum contributions from all vertices
    for (const auto& v : m_vertices)
    {
        Vector r = v - c;

        // Diagonal elements: Ixx = m(y^2 + z^2), Iyy = m(x^2 + z^2), Izz = m(x^2 + y^2)
        I(0, 0) += r[1] * r[1] + r[2] * r[2];
        I(1, 1) += r[0] * r[0] + r[2] * r[2];
        I(2, 2) += r[0] * r[0] + r[1] * r[1];

        // Off-diagonal elements (should be zero for symmetric shapes but computed for completeness)
        I(0, 1) -= r[0] * r[1];
        I(0, 2) -= r[0] * r[2];
        I(1, 2) -= r[1] * r[2];
    }

    // Average over vertices and scale by mass
    real_t n = static_cast<real_t>(m_vertices.size());
    real_t scale = m / n;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            I(i, j) *= scale;

    // Ensure symmetry
    I(1, 0) = I(0, 1);
    I(2, 0) = I(0, 2);
    I(2, 1) = I(1, 2);

    return I;
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

std::vector<Vertex> Polytope::vertices() const
{
    std::vector<Vertex> verts;
    for (const auto& v : m_vertices)
    {
        // Use normalized position as normal (points outward from center)
        Vector normal = v / norm(v);
        verts.push_back(Vertex(v, m_color, normal));
    }
    return verts;
}

std::vector<TriangleVertexIndices> Polytope::faces() const
{
    std::vector<Vertex> verts = vertices();

    if (verts.size() < 4)
        return {};

    std::vector<TriangleVertexIndices> face_indices;

    // Simplified: create triangles using first 3 vertices as base and fan from others
    for (size_t i = 3; i < verts.size(); ++i)
    {
        face_indices.push_back({0, 1, static_cast<uint_t>(i)});
        face_indices.push_back({0, static_cast<uint_t>(i), 2});
        face_indices.push_back({1, 2, static_cast<uint_t>(i)});
    }

    return face_indices;
}

Mesh Polytope::mesh() const
{
    return Mesh(vertices(), faces());
}

}  // namespace sim
