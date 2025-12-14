/*****************************************************************************************
 * @file: Mesh.hpp
 *
 * @brief: This file provides a mesh class for rendering.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include <array>
#include <vector>

#include "Camera.hpp"
#include "Color.hpp"
#include "PixelCoordinate.hpp"
#include "Transform.hpp"
#include "Vector.hpp"
#include "Vertex.hpp"
#include "types.hpp"

namespace sim
{

using TriangleVertexIndices = std::array<uint_t, 3>;

/**
 * @brief A proper Mesh class storing vertices and triangle indices.
 *
 * @details Vertices contain position (world space), color, normal, and texture coordinates.
 * Triangles are defined by indices into the vertex array.
 */
class Mesh
{
  protected:
    /**
     * @brief Vertices of the mesh.
     */
    std::vector<Vertex> m_vertices;

    /**
     * @brief Indices of the vertices belonging to each triangle.
     */
    std::vector<TriangleVertexIndices> m_triangle_indices;

  public:
    /**
     * @brief Default constructor - creates empty mesh.
     */
    Mesh();

    /**
     * @brief Construct a mesh from vertices and indices.
     * @param vertices Vector of Vertex objects.
     * @param indices Vector of triangle indices (each is 3 vertex indices).
     * @throws std::out_of_range If indices are invalid.
     */
    Mesh(const std::vector<Vertex>& vertices, const std::vector<TriangleVertexIndices>& indices);

    /**
     * @brief Get the vertices.
     */
    const std::vector<Vertex>& getVertices() const;

    /**
     * @brief Get mutable access to vertices.
     */
    std::vector<Vertex>& getVertices();

    /**
     * @brief Get the face indices.
     */
    const std::vector<TriangleVertexIndices>& getFaceIndices() const;

    /**
     * @brief Get number of vertices.
     */
    size_t getVertexCount() const;

    /**
     * @brief Get number of triangles.
     */
    size_t getTriangleCount() const;

    /**
     * @brief Set all vertex colors to a uniform color.
     */
    void setUniformColor(const Color& color);

    /**
     * @brief Add a vertex to the mesh.
     */
    void addVertex(const Vertex& vertex);

    /**
     * @brief Add a triangle to the mesh.
     * @param v0, v1, v2 The three vertices of the triangle.
     */
    void addTriangle(const Vertex& v0, const Vertex& v1, const Vertex& v2);

    /**
     * @brief Add a triangle by indices.
     * @param i0, i1, i2 Indices of the three vertices.
     */
    void addTriangleByIndices(uint_t i0, uint_t i1, uint_t i2);

    /**
     * @brief Reserve space for vertices and triangles.
     */
    void reserve(size_t vertex_count, size_t triangle_count);

    /*************************************************************************************
     * Static methods/constructors for various shapes
     ************************************************************************************/

    /**
     * @brief Create a cuboid shaped mesh.
     * @param width Width of the cuboid (in x-direction).
     * @param depth Depth of the cuboid (in y-direction).
     * @param height Height of the cuboid (in z-direction).
     * @param color Color of the cuboid.
     * @return A cuboid shaped mesh.
     */
    static Mesh Box(real_t width, real_t depth, real_t height, const Color& color);

    /**
     * @brief Create a cube shaped mesh.
     * @param side_length Side length of the cube.
     * @param color Color of the cube.
     * @return A cube shaped mesh.
     */
    static Mesh Cube(real_t side_length, const Color& color);

    /**
     * @brief Create a spherical mesh.
     * @param radius Radius of the sphere.
     * @param color Color of the sphere.
     * @param delta Angular resolution for the mesh.
     * @return A mesh approximation of a sphere.
     */
    static Mesh Sphere(real_t radius, const Color& color, real_t delta = M_PI / 16);

    /**
     * @brief Create a cylindrical mesh.
     * @param radius Radius of the cylinder.
     * @param height Height of the cylinder in z-direction.
     * @param color Color of the cylinder.
     * @param delta Angular resolution of the mesh.
     * @return A mesh approximation of a sphere.
     */
    static Mesh Cylinder(real_t radius, real_t height, const Color& color, real_t delta = M_PI / 16);

    /**
     * @brief Create a pyramid mesh with square base.
     * @param base_length Side length of the square base.
     * @param height Height of the pyramid in z-direction.
     * @param color Color of the pyramid.
     * @return A pyramid shaped mesh.
     */
    static Mesh Pyramid(real_t base_length, real_t height, const Color& color);

    /**
     * @brief Create a conical mesh.
     * @param radius Radius of the cone base.
     * @param height Height of the cone in z-direction.
     * @param color Color of the cone.
     * @param delta Angular resolution of the mesh.
     * @return A mesh approximation of a cone.
     */
    static Mesh Cone(real_t radius, real_t height, const Color& color, real_t delta = M_PI / 16);

    /*************************************************************************************
     * Computation of physical properties of a mesh
     ************************************************************************************/

    /**
     * @brief Computes the volume of a mesh.
     */
    real_t volume();

    /**
     * @brief Computes the mass of a mesh based on a uniform density.
     */
    real_t mass(real_t density);

    /**
     * @brief Computes the inertia matrix of a mesh based on a uniform density.
     */
    Matrix inertia(real_t density);

    /**
     * @brief Computes the center of mass of a mesh, assuming uniform density.
     */
    Vector center_of_mass();
};

/**
 * @brief Apply a transformation to a mesh.
 */
Mesh operator*(const Transform& transform, const Mesh& mesh);

/**
 * @brief Cull triangles that are completely outside the camera's view frustum.
 */
Mesh cullMesh(const Mesh& input_mesh, const Camera& camera);

}  // namespace sim