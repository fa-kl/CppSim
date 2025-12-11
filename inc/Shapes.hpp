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
#include "Vector.hpp"
#include "Vertex.hpp"
#include "types.hpp"

namespace sim
{

struct Color;

using Triangle = std::array<Vertex, 3>;
using Mesh = std::vector<Triangle>;

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
   * @brief Get the mesh of the sphere.
   * @return Mesh.
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

}  // namespace sim
