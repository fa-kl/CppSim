/*****************************************************************************************
 * @file: Shapes.hpp
 *
 * @brief: This file provides geometric shape classes for 2D physics simulation.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include <array>
#include <vector>

#include "Color.hpp"
#include "Mesh.hpp"
#include "Vector.hpp"
#include "Vertex.hpp"
#include "types.hpp"

namespace sim
{

struct Color;

/**
 * @brief Abstract base class for all geometric 3D shapes.
 *
 * @details All concrete shapes must implement:
 * - volume() - Compute the volume
 * - centroid() - Compute the centroid (center of mass for uniform density)
 * - inertia(ρ) - Compute moment of inertia about centroid for density ρ (3D meaning)
 * - vertices() - Get vertices for visualization and collision detection (vertices)
 * - mesh() - Get the mesh for visualizations.
 * - containsPoint(point) - Check if a point is inside the shape (3D)
 */
class Shape
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~Shape() = default;

    /**
     * @brief Compute the volume of the 3D shape.
     *
     * @return Volume of the shape.
     */
    virtual real_t volume() const = 0;

    /**
     * @brief Compute the centroid of the shape.
     *
     * @details Returns the center of mass for uniform density in shape-local coordinates.
     *
     * @return Centroid position vector.
     */
    virtual Vector centroid() const = 0;

    /**
     * @brief Compute the (scalar) moment of inertia about the centroid.
     *
     * @param rho Material density (mass per unit volume).
     * @return Scalar moment of inertia (for uniform density) about principal axes - for simple shapes this returns a
     * representative value.
     */
    virtual real_t inertia(real_t rho) const = 0;

    /**
     * @brief Get representative vertices of the shape for visualization.
     *
     * @details Returns vertices in shape-local coordinates. For curved shapes,
     * returns an approximation (e.g. tessellated sphere).
     *
     * @return Vector of 3D vertex positions.
     */
    virtual std::vector<Vertex> vertices() const = 0;

    /**
     * @brief Get mesh for visualization.
     *
     * @details Returns Mesh based on the shape's vertices.
     *
     * @return Mesh.
     */
    virtual Mesh mesh() const = 0;

    /**
     * @brief Check if a point is contained within the shape.
     *
     * @param point Point in shape-local coordinates.
     * @return True if point is inside the shape, false otherwise.
     */
    virtual bool_t containsPoint(const Vector& point) const = 0;
};

/**
 * @brief A sphere shape defined by its radius.
 *
 * @details Represents a sphere centered at the origin in local coordinates.
 */
class Sphere : public Shape
{
  protected:
    /**
     * @brief Radius of the sphere.
     */
    real_t m_radius;

    /**
     * @brief Color of the sphere.
     */
    Color m_color;

  public:
    /**
     * @brief Construct a sphere with given radius.
     *
     * @param radius Radius of the sphere (must be > 0).
     * @param color Color of the sphere.
     * @throws std::invalid_argument If radius is not positive.
     */
    explicit Sphere(real_t radius, Color color);

    /**
     * @brief Get the radius of the sphere.
     *
     * @return Radius.
     */
    real_t radius() const;

    /**
     * @brief Compute the volume of the sphere.
     *
     * @details Volume = (4/3) π r^3
     *
     * @return Volume of the sphere.
     */
    real_t volume() const override;

    /**
     * @brief Compute the centroid of the sphere.
     *
     * @details Always at origin for shape-local coordinates.
     *
     * @return Centroid position vector (0, 0, 0).
     */
    Vector centroid() const override;

    /**
     * @brief Compute the moment of inertia about the centroid.
     *
     * @details For a solid sphere: I = (2/5) * m * r^2 where m = ρ * V
     *
     * @param rho Material density (mass per unit volume).
     * @return Scalar representative moment of inertia.
     */
    real_t inertia(real_t rho) const override;

    /**
     * @brief Get vertices approximating the sphere.
     *
     * @details Returns a tessellation of the sphere surface using angular step `delta`.
     *
     * @param delta Angular step size in radians (default: π/8).
     * @return Vertices approximating the sphere.
     */
    std::vector<Vertex> vertices(real_t delta) const;

    /**
     * @brief Get vertices approximating the sphere (override for Shape interface).
     *
     * @return Vertices approximating the sphere.
     */
    std::vector<Vertex> vertices() const override;

    /**
     * @brief Get the mesh of the sphere with custom tessellation.
     * @param delta Angular step size in radians.
     * @return Mesh.
     */
    Mesh mesh(real_t delta) const;

    /**
     * @brief Get the mesh of the sphere (override for Shape interface).
     * @return Mesh with default tessellation.
     */
    Mesh mesh() const override;

    /**
     * @brief Check if a point is contained within the sphere.
     *
     * @param point Point in shape-local coordinates.
     * @return True if point is inside the sphere, false otherwise.
     */
    bool_t containsPoint(const Vector& point) const override;
};

/**
 * @brief A cuboid shape defined by width, height and depth.
 *
 * @details Represents an axis-aligned cuboid centered at the origin in local coordinates.
 */
class Box : public Shape
{
  protected:
    /**
     * @brief Width of the cuboid.
     */
    real_t m_width;

    /**
     * @brief Height of the cuboid.
     */
    real_t m_height;

    /**
     * @brief Depth of the cuboid.
     */
    real_t m_depth;

    /**
     * @brief Color of the cuboid.
     */
    Color m_color;

  public:
    /**
     * @brief Construct a cuboid with given width, height and depth.
     *
     * @param width Width of the box (must be > 0).
     * @param height Height of the box (must be > 0).
     * @param depth Depth of the box (must be > 0).
     * @param color Color of the box.
     * @throws std::invalid_argument If any dimension is not positive.
     */
    Box(real_t width, real_t height, real_t depth, Color color);

    /**
     * @brief Construct a cube with equal side lengths.
     *
     * @param size Side length of the cube (must be > 0).
     * @param color Color of the box.
     * @throws std::invalid_argument If size is not positive.
     */
    explicit Box(real_t size, Color color);

    real_t width() const;
    real_t height() const;
    real_t depth() const;

    /**
     * @brief Compute the volume of the cuboid.
     *
     * @return Volume = w * h * d.
     */
    real_t volume() const override;

    Vector centroid() const override;

    /**
     * @brief Compute the (representative) moment of inertia about the centroid.
     *
     * @details Principal inertias for a cuboid are:
     * I_x = (1/12) m (h^2 + d^2)
     * I_y = (1/12) m (w^2 + d^2)
     * I_z = (1/12) m (w^2 + h^2)
     * This function returns the average of the three principal values as a scalar representative.
     *
     * @param rho Material density (mass per unit volume).
     * @return Representative scalar inertia.
     */
    real_t inertia(real_t rho) const override;

    /**
     * @brief Get the eight corner vertices of the cuboid.
     *
     * @return Vertices of the box.
     */
    std::vector<Vertex> vertices() const override;

    /**
     * @brief Get the mesh of the cuboid.
     * @return Mesh.
     */
    Mesh mesh() const override;

    bool_t containsPoint(const Vector& point) const override;
};

/**
 * @brief A cylinder shape defined by radius and height.
 *
 * @details Represents a cylinder centered at the origin with axis along the z-axis.
 */
class Cylinder : public Shape
{
  protected:
    /**
     * @brief Radius of the cylinder.
     */
    real_t m_radius;

    /**
     * @brief Height of the cylinder.
     */
    real_t m_height;

    /**
     * @brief Color of the cylinder.
     */
    Color m_color;

  public:
    /**
     * @brief Construct a cylinder with given radius and height.
     *
     * @param radius Radius of the cylinder (must be > 0).
     * @param height Height of the cylinder (must be > 0).
     * @param color Color of the cylinder.
     * @throws std::invalid_argument If radius or height is not positive.
     */
    Cylinder(real_t radius, real_t height, Color color);

    /**
     * @brief Get the radius of the cylinder.
     * @return Radius.
     */
    real_t radius() const;

    /**
     * @brief Get the height of the cylinder.
     * @return Height.
     */
    real_t height() const;

    /**
     * @brief Compute the volume of the cylinder.
     * @details Volume = π r^2 h
     * @return Volume of the cylinder.
     */
    real_t volume() const override;

    Vector centroid() const override;

    /**
     * @brief Compute the moment of inertia about the centroid.
     * @details For a solid cylinder: I_z = (1/2) m r^2, I_x = I_y = (1/12) m (3r^2 + h^2)
     * @param rho Material density (mass per unit volume).
     * @return Representative scalar inertia.
     */
    real_t inertia(real_t rho) const override;

    /**
     * @brief Get vertices approximating the cylinder.
     * @param delta Angular step size in radians.
     * @return Vertices approximating the cylinder.
     */
    std::vector<Vertex> vertices(real_t delta) const;

    std::vector<Vertex> vertices() const override;

    /**
     * @brief Get the mesh of the cylinder with custom tessellation.
     * @param delta Angular step size in radians.
     * @return Mesh.
     */
    Mesh mesh(real_t delta) const;

    Mesh mesh() const override;

    bool_t containsPoint(const Vector& point) const override;
};

/**
 * @brief A pyramid shape with a square base.
 *
 * @details Represents a pyramid centered at the origin with base in the xy-plane.
 */
class Pyramid : public Shape
{
  protected:
    /**
     * @brief Base width/length of the pyramid.
     */
    real_t m_base;

    /**
     * @brief Height of the pyramid.
     */
    real_t m_height;

    /**
     * @brief Color of the pyramid.
     */
    Color m_color;

  public:
    /**
     * @brief Construct a pyramid with given base size and height.
     *
     * @param base Base width/length (must be > 0).
     * @param height Height of the pyramid (must be > 0).
     * @param color Color of the pyramid.
     * @throws std::invalid_argument If base or height is not positive.
     */
    Pyramid(real_t base, real_t height, Color color);

    /**
     * @brief Get the base size of the pyramid.
     * @return Base width/length.
     */
    real_t base() const;

    /**
     * @brief Get the height of the pyramid.
     * @return Height.
     */
    real_t height() const;

    /**
     * @brief Compute the volume of the pyramid.
     * @details Volume = (1/3) * base^2 * height
     * @return Volume of the pyramid.
     */
    real_t volume() const override;

    Vector centroid() const override;

    /**
     * @brief Compute the moment of inertia about the centroid.
     * @param rho Material density (mass per unit volume).
     * @return Representative scalar inertia.
     */
    real_t inertia(real_t rho) const override;

    std::vector<Vertex> vertices() const override;

    Mesh mesh() const override;

    bool_t containsPoint(const Vector& point) const override;
};

/**
 * @brief A cone shape defined by radius and height.
 *
 * @details Represents a cone centered at the origin with base in the xy-plane and apex along +z.
 */
class Cone : public Shape
{
  protected:
    /**
     * @brief Radius of the cone base.
     */
    real_t m_radius;

    /**
     * @brief Height of the cone.
     */
    real_t m_height;

    /**
     * @brief Color of the cone.
     */
    Color m_color;

  public:
    /**
     * @brief Construct a cone with given radius and height.
     *
     * @param radius Radius of the base (must be > 0).
     * @param height Height of the cone (must be > 0).
     * @param color Color of the cone.
     * @throws std::invalid_argument If radius or height is not positive.
     */
    Cone(real_t radius, real_t height, Color color);

    /**
     * @brief Get the radius of the cone.
     * @return Radius.
     */
    real_t radius() const;

    /**
     * @brief Get the height of the cone.
     * @return Height.
     */
    real_t height() const;

    /**
     * @brief Compute the volume of the cone.
     * @details Volume = (1/3) π r^2 h
     * @return Volume of the cone.
     */
    real_t volume() const override;

    Vector centroid() const override;

    /**
     * @brief Compute the moment of inertia about the centroid.
     * @param rho Material density (mass per unit volume).
     * @return Representative scalar inertia.
     */
    real_t inertia(real_t rho) const override;

    /**
     * @brief Get vertices approximating the cone.
     * @param delta Angular step size in radians.
     * @return Vertices approximating the cone.
     */
    std::vector<Vertex> vertices(real_t delta) const;

    std::vector<Vertex> vertices() const override;

    /**
     * @brief Get the mesh of the cone with custom tessellation.
     * @param delta Angular step size in radians.
     * @return Mesh.
     */
    Mesh mesh(real_t delta) const;

    Mesh mesh() const override;

    bool_t containsPoint(const Vector& point) const override;
};

/**
 * @brief A convex polytope defined by arbitrary vertices.
 *
 * @details Represents a convex polyhedron defined by its vertices.
 * The convex hull of the vertices is computed and triangulated.
 */
class Polytope : public Shape
{
  protected:
    /**
     * @brief Vertices defining the polytope.
     */
    std::vector<Vector> m_vertices;

    /**
     * @brief Color of the polytope.
     */
    Color m_color;

    /**
     * @brief Precomputed centroid.
     */
    mutable Vector m_centroid;

    /**
     * @brief Whether centroid has been computed.
     */
    mutable bool m_centroid_computed;

  public:
    /**
     * @brief Construct a polytope from a set of vertices.
     *
     * @param vertices Vector of 3D vertices defining the polytope.
     * @param color Color of the polytope.
     * @throws std::invalid_argument If fewer than 4 vertices provided.
     */
    Polytope(const std::vector<Vector>& vertices, Color color);

    /**
     * @brief Get the vertices of the polytope.
     * @return Vector of vertex positions.
     */
    const std::vector<Vector>& getVertices() const;

    /**
     * @brief Compute the volume of the polytope.
     * @details Uses divergence theorem with centroid as reference point.
     * @return Volume of the polytope.
     */
    real_t volume() const override;

    Vector centroid() const override;

    /**
     * @brief Compute the moment of inertia about the centroid.
     * @param rho Material density (mass per unit volume).
     * @return Representative scalar inertia.
     */
    real_t inertia(real_t rho) const override;

    std::vector<Vertex> vertices() const override;

    Mesh mesh() const override;

    bool_t containsPoint(const Vector& point) const override;
};

}  // namespace sim
