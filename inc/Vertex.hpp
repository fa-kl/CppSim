/*****************************************************************************************
 * @file: Vertex.hpp
 *
 * @brief: Vertex struct for rendering with position, color, normal, and texture coordinates.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include "Color.hpp"
#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

/**
 * @brief A vertex struct for 3D rendering.
 *
 * @details Stores position (world space), color, normal vector, and texture coordinates.
 * This is designed to be compatible with OBJ file format and modern rendering pipelines.
 */
struct Vertex
{
    /**
     * @brief 3D position in world space coordinates.
     */
    Vector position;

    /**
     * @brief Vertex color (RGBA).
     */
    Color color;

    /**
     * @brief Normal vector (3D, should be normalized).
     */
    Vector normal;

    /**
     * @brief Construct a vertex with position, color, and normal.
     * @param p 3D position vector.
     * @param c Vertex color.
     * @param n 3D normal vector.
     */
    Vertex(const Vector& p, const Color& c, const Vector& n);
};

}  // namespace sim