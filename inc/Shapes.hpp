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

#include <vector>

#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

/**
 * @brief Abstract base class for all geometric shapes.
 *
 * @details All concrete shapes must implement:
 * - area() - Compute the area
 * - centroid() - Compute the centroid (center of mass for uniform density)
 * - inertia(ρ) - Compute moment of inertia about centroid for density ρ
 * - vertices() - Get vertices for visualization and collision detection
 * - containsPoint(point) - Check if a point is inside the shape
 */
class Shape
{
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~Shape() = default;

  /**
   * @brief Compute the area of the shape.
   *
   * @return Area of the shape.
   */
  virtual real_t area() const = 0;

  /**
   * @brief Compute the centroid of the shape.
   *
   * @details Returns the center of mass for uniform density in shape-local coordinates.
   *
   * @return Centroid position vector.
   */
  virtual Vector centroid() const = 0;

  /**
   * @brief Compute the moment of inertia about the centroid.
   *
   * @param rho Material density (mass per unit area).
   * @return Moment of inertia about the z-axis through centroid.
   */
  virtual real_t inertia(real_t rho) const = 0;

  /**
   * @brief Get the vertices of the shape for visualization.
   *
   * @details Returns vertices in counter-clockwise order in shape-local coordinates.
   * For curved shapes, returns an approximation.
   *
   * @return Vector of vertex positions.
   */
  virtual std::vector<Vector> vertices() const = 0;

  /**
   * @brief Check if a point is contained within the shape.
   *
   * @param point Point in shape-local coordinates.
   * @return True if point is inside the shape, false otherwise.
   */
  virtual bool_t containsPoint(const Vector& point) const = 0;
};

/**
 * @brief A circle shape defined by its radius.
 *
 * @details Represents a circular shape centered at the origin in local coordinates.
 */
class Circle : public Shape
{
protected:
  /**
   * @brief Radius of the circle.
   */
  real_t m_radius;

public:
  /**
   * @brief Construct a circle with given radius.
   *
   * @param radius Radius of the circle (must be > 0).
   * @throws std::invalid_argument If radius is not positive.
   */
  explicit Circle(real_t radius);

  /**
   * @brief Get the radius of the circle.
   *
   * @return Radius.
   */
  real_t radius() const;

  /**
   * @brief Compute the area of the circle.
   *
   * @details Area = πr²
   *
   * @return Area of the circle.
   */
  real_t area() const override;

  /**
   * @brief Compute the centroid of the circle.
   *
   * @details Always at origin for shape-local coordinates.
   *
   * @return Centroid position vector (0, 0).
   */
  Vector centroid() const override;

  /**
   * @brief Compute the moment of inertia about the centroid.
   *
   * @details For a circle: I = (1/2) * m * r² where m = ρ * A
   *
   * @param rho Material density (mass per unit area).
   * @return Moment of inertia about the z-axis through centroid.
   */
  real_t inertia(real_t rho) const override;

  /**
   * @brief Get vertices approximating the circle.
   *
   * @details Returns vertices in counter-clockwise order starting from the rightmost point.
   *
   * @param delta_phi Angular step size in radians (default: π/8).
   * @return Vector of vertex positions approximating the circle.
   */
  std::vector<Vector> vertices(real_t delta_phi) const;

  /**
   * @brief Get vertices approximating the circle (override for Shape interface).
   *
   * @return Vector of vertex positions approximating the circle.
   */
  std::vector<Vector> vertices() const override;

  /**
   * @brief Check if a point is contained within the circle.
   *
   * @param point Point in shape-local coordinates.
   * @return True if point is inside the circle, false otherwise.
   */
  bool_t containsPoint(const Vector& point) const override;
};

/**
 * @brief A rectangular box shape defined by width and height.
 *
 * @details Represents a rectangle centered at the origin in local coordinates,
 * with sides parallel to the coordinate axes.
 */
class Box : public Shape
{
protected:
  /**
   * @brief Width of the box (along x-axis).
   */
  real_t m_width;

  /**
   * @brief Height of the box (along y-axis).
   */
  real_t m_height;

public:
  /**
   * @brief Construct a box with given width and height.
   *
   * @param width Width of the box (must be > 0).
   * @param height Height of the box (must be > 0).
   * @throws std::invalid_argument If width or height is not positive.
   */
  Box(real_t width, real_t height);

  /**
   * @brief Construct a square box with equal width and height.
   *
   * @param size Side length of the square (must be > 0).
   * @throws std::invalid_argument If size is not positive.
   */
  explicit Box(real_t size);

  /**
   * @brief Get the width of the box.
   *
   * @return Width.
   */
  real_t width() const;

  /**
   * @brief Get the height of the box.
   *
   * @return Height.
   */
  real_t height() const;

  /**
   * @brief Compute the area of the box.
   *
   * @details Area = width * height
   *
   * @return Area of the box.
   */
  real_t area() const override;

  /**
   * @brief Compute the centroid of the box.
   *
   * @details Always at origin for shape-local coordinates.
   *
   * @return Centroid position vector (0, 0).
   */
  Vector centroid() const override;

  /**
   * @brief Compute the moment of inertia about the centroid.
   *
   * @details For a rectangle: I = (1/12) * m * (w² + h²) where m = ρ * A
   *
   * @param rho Material density (mass per unit area).
   * @return Moment of inertia about the z-axis through centroid.
   */
  real_t inertia(real_t rho) const override;

  /**
   * @brief Get the four corner vertices of the box.
   *
   * @details Returns vertices in counter-clockwise order starting from bottom-left.
   *
   * @return Vector of vertex positions.
   */
  std::vector<Vector> vertices() const override;

  /**
   * @brief Check if a point is contained within the box.
   *
   * @param point Point in shape-local coordinates.
   * @return True if point is inside the box, false otherwise.
   */
  bool_t containsPoint(const Vector& point) const override;
};

/**
 * @brief A convex polygon shape defined by vertices.
 *
 * @details Represents a convex polygon with vertices in counter-clockwise order.
 * The polygon is assumed to be convex; centroid and inertia calculations are only
 * valid for convex polygons.
 */
class Polygon : public Shape
{
protected:
  /**
   * @brief Vertices in counter-clockwise order.
   */
  std::vector<Vector> m_vertices;

public:
  /**
   * @brief Construct a polygon from vertices.
   *
   * @param vertices Vertices in counter-clockwise order (must have ≥ 3 vertices).
   * @throws std::invalid_argument If less than 3 vertices provided.
   */
  explicit Polygon(const std::vector<Vector>& vertices);

  /**
   * @brief Construct a polygon from an initializer list.
   *
   * @param vertices Initializer list of vertices.
   */
  Polygon(std::initializer_list<Vector> vertices);

  /**
   * @brief Get the vertices of the polygon.
   *
   * @return Vector of vertex positions.
   */
  const std::vector<Vector>& getVertices() const;

  /**
   * @brief Compute the area of the polygon using the shoelace formula.
   *
   * @details Area = (1/2) * |∑(xᵢ * yᵢ₊₁ - xᵢ₊₁ * yᵢ)|
   * Returns positive area for counter-clockwise vertex ordering.
   *
   * @return Area of the polygon.
   */
  real_t area() const override;

  /**
   * @brief Compute the centroid of the polygon.
   *
   * @details For a polygon with vertices (xᵢ, yᵢ):
   * - Cₓ = (1/6A) * ∑(xᵢ + xᵢ₊₁)(xᵢyᵢ₊₁ - xᵢ₊₁yᵢ)
   * - Cᵧ = (1/6A) * ∑(yᵢ + yᵢ₊₁)(xᵢyᵢ₊₁ - xᵢ₊₁yᵢ)
   *
   * @return Centroid position vector.
   */
  Vector centroid() const override;

  /**
   * @brief Compute the moment of inertia about the centroid.
   *
   * @details Uses the formula:
   * I = ρ * ∑[(xᵢyᵢ₊₁ - xᵢ₊₁yᵢ) * (xᵢ² + xᵢxᵢ₊₁ + xᵢ₊₁² + yᵢ² + yᵢyᵢ₊₁ + yᵢ₊₁²)] / 12
   * Then applies parallel axis theorem to shift to centroid.
   *
   * @param rho Material density (mass per unit area).
   * @return Moment of inertia about the z-axis through centroid.
   */
  real_t inertia(real_t rho) const override;

  /**
   * @brief Get the vertices of the polygon.
   *
   * @details Returns vertices in counter-clockwise order.
   *
   * @return Vector of vertex positions.
   */
  std::vector<Vector> vertices() const override;

  /**
   * @brief Check if a point is contained within the polygon.
   *
   * @details Uses the winding number algorithm for point-in-polygon test.
   *
   * @param point Point in shape-local coordinates.
   * @return True if point is inside the polygon, false otherwise.
   */
  bool_t containsPoint(const Vector& point) const override;
};

}  // namespace sim
