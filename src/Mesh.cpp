/*****************************************************************************************
 * @file: Mesh.cpp
 *
 * @brief: This file implements a mesh class for rendering.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Mesh.hpp"

#include <cmath>
#include <stdexcept>

namespace sim
{

Mesh::Mesh() : m_vertices(), m_triangle_indices() {}

Mesh::Mesh(const std::vector<Vertex>& vertices, const std::vector<TriangleVertexIndices>& indices)
    : m_vertices(vertices), m_triangle_indices(indices)
{
    for (const TriangleVertexIndices& face : m_triangle_indices)
    {
        for (size_t idx : face)
        {
            if (idx >= m_vertices.size())
            {
                throw std::out_of_range("Face index out of range");
            }
        }
    }
}

const std::vector<Vertex>& Mesh::getVertices() const
{
    return m_vertices;
}

std::vector<Vertex>& Mesh::getVertices()
{
    return m_vertices;
}

const std::vector<TriangleVertexIndices>& Mesh::getFaceIndices() const
{
    return m_triangle_indices;
}

size_t Mesh::getVertexCount() const
{
    return m_vertices.size();
}

size_t Mesh::getTriangleCount() const
{
    return m_triangle_indices.size();
}

void Mesh::setUniformColor(const Color& color)
{
    for (auto& vertex : m_vertices)
    {
        vertex.color = color;
    }
}

void Mesh::addVertex(const Vertex& vertex)
{
    m_vertices.push_back(vertex);
}

void Mesh::addTriangle(const Vertex& v0, const Vertex& v1, const Vertex& v2)
{
    uint_t idx0 = m_vertices.size();
    uint_t idx1 = idx0 + 1;
    uint_t idx2 = idx0 + 2;

    m_vertices.push_back(v0);
    m_vertices.push_back(v1);
    m_vertices.push_back(v2);

    m_triangle_indices.push_back({idx0, idx1, idx2});
}

void Mesh::addTriangleByIndices(uint_t i0, uint_t i1, uint_t i2)
{
    if (i0 >= m_vertices.size() || i1 >= m_vertices.size() || i2 >= m_vertices.size())
    {
        throw std::out_of_range("Triangle index out of range");
    }
    m_triangle_indices.push_back({i0, i1, i2});
}

void Mesh::reserve(size_t vertex_count, size_t triangle_count)
{
    m_vertices.reserve(vertex_count);
    m_triangle_indices.reserve(triangle_count);
}

/*================================================================================================
 * Static methods for creating various shapes
 ================================================================================================*/

Mesh Mesh::Box(real_t width, real_t depth, real_t height, const Color& color)
{
    Mesh mesh;
    mesh.reserve(8, 12);
    real_t hw = width / 2.0;
    real_t hd = depth / 2.0;
    real_t hh = height / 2.0;
    Vector corners[8] = {
        Vector{-hw, -hd, -hh},  // 1: left-front-bottom
        Vector{hw, -hd, -hh},   // 2: right-front-bottom
        Vector{hw, hd, -hh},    // 3: right-back-bottom
        Vector{-hw, hd, -hh},   // 4: left-back-bottom
        Vector{-hw, -hd, hh},   // 5: left-front-top
        Vector{hw, -hd, hh},    // 6: right-front-top
        Vector{hw, hd, hh},     // 7: right-back-top
        Vector{-hw, hd, hh}     // 8: left-back-top
    };
    for (size_t i = 0; i < 8; ++i)
    {
        mesh.addVertex(Vertex(corners[i], color, normalize(corners[i])));
    }
    // Bottom face
    mesh.addTriangleByIndices(0, 2, 1);
    mesh.addTriangleByIndices(0, 3, 2);
    // Top face
    mesh.addTriangleByIndices(4, 5, 6);
    mesh.addTriangleByIndices(4, 6, 7);
    // Front face
    mesh.addTriangleByIndices(0, 1, 5);
    mesh.addTriangleByIndices(0, 5, 4);
    // Back face
    mesh.addTriangleByIndices(2, 3, 7);
    mesh.addTriangleByIndices(2, 7, 6);
    // Left face
    mesh.addTriangleByIndices(3, 0, 4);
    mesh.addTriangleByIndices(3, 4, 7);
    // Right face
    mesh.addTriangleByIndices(1, 2, 6);
    mesh.addTriangleByIndices(1, 6, 5);

    return mesh;
}

Mesh Mesh::Cube(real_t side_length, const Color& color)
{
    return Box(side_length, side_length, side_length, color);
}

Mesh Mesh::Sphere(real_t radius, const Color& color, real_t delta)
{
    Mesh mesh;

    // Create sphere using spherical coordinates
    // theta: azimuthal angle (0 to 2*pi)
    // phi: polar angle (0 to pi)

    std::vector<std::vector<uint_t>> grid_indices;

    // Generate vertices
    for (real_t phi = 0; phi <= M_PI; phi += delta)
    {
        std::vector<uint_t> row_indices;
        for (real_t theta = 0; theta < 2 * M_PI; theta += delta)
        {
            real_t x = radius * std::sin(phi) * std::cos(theta);
            real_t y = radius * std::sin(phi) * std::sin(theta);
            real_t z = radius * std::cos(phi);

            Vector position{x, y, z};
            Vector normal = normalize(position);  // normalized position

            mesh.addVertex(Vertex(position, color, normal));
            row_indices.push_back(mesh.getVertexCount() - 1);
        }
        grid_indices.push_back(row_indices);
    }

    // Create triangles from grid
    for (size_t i = 0; i < grid_indices.size() - 1; ++i)
    {
        for (size_t j = 0; j < grid_indices[i].size(); ++j)
        {
            size_t next_j = (j + 1) % grid_indices[i].size();

            uint_t v0 = grid_indices[i][j];
            uint_t v1 = grid_indices[i][next_j];
            uint_t v2 = grid_indices[i + 1][next_j];
            uint_t v3 = grid_indices[i + 1][j];

            mesh.addTriangleByIndices(v0, v1, v2);
            mesh.addTriangleByIndices(v0, v2, v3);
        }
    }

    return mesh;
}

Mesh Mesh::Cylinder(real_t radius, real_t height, const Color& color, real_t delta)
{
    Mesh mesh;

    real_t half_height = height / 2.0;
    uint_t num_segments = static_cast<uint_t>(2 * M_PI / delta);

    // Center vertices for caps
    uint_t bottom_center_idx = 0;
    mesh.addVertex(Vertex(Vector{0, 0, -half_height}, color, Vector{0, 0, -1}));

    uint_t top_center_idx = 1;
    mesh.addVertex(Vertex(Vector{0, 0, half_height}, color, Vector{0, 0, 1}));

    // Generate circle vertices (shared between caps and sides)
    // We'll use averaged normals for smooth shading
    for (size_t i = 0; i < num_segments; ++i)
    {
        real_t theta = static_cast<real_t>(i) * delta;
        real_t x = radius * std::cos(theta);
        real_t y = radius * std::sin(theta);

        // Radial normal for smooth shading on sides
        Vector radial_normal{x / radius, y / radius, 0};

        // Bottom circle vertex
        mesh.addVertex(Vertex(Vector{x, y, -half_height}, color, radial_normal));
        // Top circle vertex
        mesh.addVertex(Vertex(Vector{x, y, half_height}, color, radial_normal));
    }

    // Create bottom cap triangles
    for (size_t i = 0; i < num_segments; ++i)
    {
        uint_t current_bottom = 2 + 2 * i;
        uint_t next_bottom = 2 + 2 * ((i + 1) % num_segments);
        mesh.addTriangleByIndices(bottom_center_idx, next_bottom, current_bottom);
    }

    // Create top cap triangles
    for (size_t i = 0; i < num_segments; ++i)
    {
        uint_t current_top = 2 + 2 * i + 1;
        uint_t next_top = 2 + 2 * ((i + 1) % num_segments) + 1;
        mesh.addTriangleByIndices(top_center_idx, current_top, next_top);
    }

    // Create side triangles
    for (size_t i = 0; i < num_segments; ++i)
    {
        uint_t current_bottom = 2 + 2 * i;
        uint_t current_top = 2 + 2 * i + 1;
        uint_t next_bottom = 2 + 2 * ((i + 1) % num_segments);
        uint_t next_top = 2 + 2 * ((i + 1) % num_segments) + 1;

        mesh.addTriangleByIndices(current_bottom, current_top, next_top);
        mesh.addTriangleByIndices(current_bottom, next_top, next_bottom);
    }

    return mesh;
}

Mesh Mesh::Pyramid(real_t base_length, real_t height, const Color& color)
{
    Mesh mesh;
    mesh.reserve(5, 6);  // 5 vertices (4 base + 1 apex), 6 triangles

    real_t half_base = base_length / 2.0;

    // Base vertices (z = 0) with normals pointing downward/outward
    mesh.addVertex(Vertex(Vector{-half_base, -half_base, 0}, color, normalize(Vector{-1, -1, -1})));
    mesh.addVertex(Vertex(Vector{half_base, -half_base, 0}, color, normalize(Vector{1, -1, -1})));
    mesh.addVertex(Vertex(Vector{half_base, half_base, 0}, color, normalize(Vector{1, 1, -1})));
    mesh.addVertex(Vertex(Vector{-half_base, half_base, 0}, color, normalize(Vector{-1, 1, -1})));

    // Apex with normal pointing upward
    mesh.addVertex(Vertex(Vector{0, 0, height}, color, Vector{0, 0, 1}));

    // Base triangles
    mesh.addTriangleByIndices(0, 2, 1);
    mesh.addTriangleByIndices(0, 3, 2);

    // Side triangles (apex is vertex 4)
    mesh.addTriangleByIndices(0, 1, 4);  // Front face
    mesh.addTriangleByIndices(1, 2, 4);  // Right face
    mesh.addTriangleByIndices(2, 3, 4);  // Back face
    mesh.addTriangleByIndices(3, 0, 4);  // Left face

    return mesh;
}

Mesh Mesh::Cone(real_t radius, real_t height, const Color& color, real_t delta)
{
    Mesh mesh;

    uint_t num_segments = static_cast<uint_t>(2 * M_PI / delta);

    // Base center
    uint_t base_center_idx = 0;
    mesh.addVertex(Vertex(Vector{0, 0, 0}, color, Vector{0, 0, -1}));

    // Apex (vertex 1)
    uint_t apex_idx = 1;
    mesh.addVertex(Vertex(Vector{0, 0, height}, color, Vector{0, 0, 1}));

    // Generate base circle vertices
    // Calculate cone side normal angle (slant)
    real_t slant = atan2(radius, height);
    for (size_t i = 0; i < num_segments; ++i)
    {
        real_t theta = static_cast<real_t>(i) * delta;
        real_t x = radius * std::cos(theta);
        real_t y = radius * std::sin(theta);

        // Normal pointing outward along the cone surface
        Vector radial{x / radius, y / radius, 0};
        Vector side_normal =
            normalize(Vector{radial[0] * std::cos(slant), radial[1] * std::cos(slant), std::sin(slant)});

        mesh.addVertex(Vertex(Vector{x, y, 0}, color, side_normal));
    }

    // Create base triangles
    for (size_t i = 0; i < num_segments; ++i)
    {
        uint_t current = 2 + i;
        uint_t next = 2 + ((i + 1) % num_segments);
        mesh.addTriangleByIndices(base_center_idx, next, current);
    }

    // Create side triangles
    for (size_t i = 0; i < num_segments; ++i)
    {
        uint_t current = 2 + i;
        uint_t next = 2 + ((i + 1) % num_segments);
        mesh.addTriangleByIndices(current, next, apex_idx);
    }

    return mesh;
}

/*================================================================================================
 * Physical property computations
 ================================================================================================*/

real_t Mesh::volume()
{
    // Use signed volume of tetrahedra from origin to each triangle
    // Sum over all triangles: V = (1/6) * sum(v0 · (v1 × v2))
    real_t total_volume = 0.0;

    for (const TriangleVertexIndices& face : m_triangle_indices)
    {
        const Vector& v0 = m_vertices[face[0]].position;
        const Vector& v1 = m_vertices[face[1]].position;
        const Vector& v2 = m_vertices[face[2]].position;

        // Signed volume of tetrahedron formed by origin and triangle
        real_t signed_vol = dot(v0, cross(v1, v2));
        total_volume += signed_vol;
    }

    return std::abs(total_volume) / 6.0;
}

real_t Mesh::mass(real_t density)
{
    return density * volume();
}

Vector Mesh::center_of_mass()
{
    // For uniform density, center of mass weighted by tetrahedron volumes
    Vector com{0, 0, 0};
    real_t total_signed_volume = 0.0;

    for (const TriangleVertexIndices& face : m_triangle_indices)
    {
        const Vector& v0 = m_vertices[face[0]].position;
        const Vector& v1 = m_vertices[face[1]].position;
        const Vector& v2 = m_vertices[face[2]].position;

        // Signed volume of tetrahedron
        real_t signed_vol = dot(v0, cross(v1, v2)) / 6.0;

        // Centroid of tetrahedron is at (v0 + v1 + v2) / 4
        Vector centroid = (v0 + v1 + v2) * 0.25;

        com = com + centroid * signed_vol;
        total_signed_volume += signed_vol;
    }

    if (std::abs(total_signed_volume) > 1e-10)
    {
        com = com * (1.0 / total_signed_volume);
    }

    return com;
}

Matrix Mesh::inertia(real_t density)
{
    // Compute inertia tensor for uniform density mesh
    // Using canonical tetrahedron inertia formulas
    Matrix I(3, 3);

    for (const TriangleVertexIndices& face : m_triangle_indices)
    {
        const Vector& v0 = m_vertices[face[0]].position;
        const Vector& v1 = m_vertices[face[1]].position;
        const Vector& v2 = m_vertices[face[2]].position;

        // Signed volume
        real_t signed_vol = dot(v0, cross(v1, v2)) / 6.0;

        // For a tetrahedron with one vertex at origin and three at v0, v1, v2
        // we compute the inertia tensor contribution
        // mass of this tetrahedron
        real_t tet_mass = density * std::abs(signed_vol);

        // Create matrix from vertices for inertia calculation
        // Using formula for tetrahedron inertia tensor
        // Note: operator() uses 1-based indexing
        for (size_t i = 1; i <= 3; ++i)
        {
            for (size_t j = 1; j <= 3; ++j)
            {
                real_t sum = 0.0;

                // Sum of products for tetrahedron vertices (0-based for Vector access)
                sum += v0[i - 1] * v0[j - 1];
                sum += v1[i - 1] * v1[j - 1];
                sum += v2[i - 1] * v2[j - 1];
                sum += v0[i - 1] * v1[j - 1] + v1[i - 1] * v0[j - 1];
                sum += v0[i - 1] * v2[j - 1] + v2[i - 1] * v0[j - 1];
                sum += v1[i - 1] * v2[j - 1] + v2[i - 1] * v1[j - 1];

                if (i == j)
                {
                    // Diagonal: I_ii = mass * (sum of squared distances - sum_ii)
                    real_t trace_contrib = 0.0;
                    for (size_t k = 1; k <= 3; ++k)
                    {
                        if (k != i)
                        {
                            real_t sum_kk = 0.0;
                            sum_kk += v0[k - 1] * v0[k - 1];
                            sum_kk += v1[k - 1] * v1[k - 1];
                            sum_kk += v2[k - 1] * v2[k - 1];
                            sum_kk += v0[k - 1] * v1[k - 1];
                            sum_kk += v0[k - 1] * v2[k - 1];
                            sum_kk += v1[k - 1] * v2[k - 1];
                            trace_contrib += sum_kk;
                        }
                    }
                    I(static_cast<index_t>(i), static_cast<index_t>(j)) = I(static_cast<index_t>(i), static_cast<index_t>(j)) + (tet_mass / 10.0) * trace_contrib;
                }
                else
                {
                    // Off-diagonal: I_ij = -mass * sum_ij
                    I(static_cast<index_t>(i), static_cast<index_t>(j)) = I(static_cast<index_t>(i), static_cast<index_t>(j)) - (tet_mass / 20.0) * sum;
                }
            }
        }
    }

    return I;
}

Mesh operator*(const Transform& transform, const Mesh& mesh)
{
    std::vector<Vertex> transformed_vertices;
    transformed_vertices.reserve(mesh.getVertexCount());

    for (const Vertex& vertex : mesh.getVertices())
    {
        Vector transformed_pos = transform * vertex.position;
        Vector transformed_normal = transform.rotation() * vertex.normal;
        Vertex transformed_vertex(transformed_pos, vertex.color, transformed_normal);
        transformed_vertices.push_back(transformed_vertex);
    }
    return Mesh(transformed_vertices, mesh.getFaceIndices());
}

Mesh cullMesh(const Mesh& input_mesh, const Camera& camera)
{
    Mesh result;
    result.reserve(input_mesh.getVertexCount(), input_mesh.getTriangleCount());
    const std::vector<Vertex>& vertices = input_mesh.getVertices();
    for (const TriangleVertexIndices& face : input_mesh.getFaceIndices())
    {
        bool any_in_view = false;
        for (size_t idx : face)
        {
            if (camera.isInView(vertices[idx].position))
            {
                any_in_view = true;
                break;
            }
        }
        if (any_in_view)
        {
            result.addTriangle(vertices[face[0]], vertices[face[1]], vertices[face[2]]);
        }
    }
    return result;
}

}  // namespace sim
