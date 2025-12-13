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
#include "Transform.hpp"
#include "Vertex.hpp"
#include "types.hpp"

namespace sim
{

struct Color;

/**
 * @brief A shortcut for triangles (3 vertices).
 */
using Triangle = std::array<Vertex, 3>;

/**
 * @brief A type for meshes.
 */
using Mesh = std::vector<Triangle>;

/**
 * @brief Apply a transformation to a mesh.
 */
Mesh operator*(const Transform& transform, const Mesh& mesh);

/**
 * @brief Cull triangles that are completely outside the camera's view frustum.
 */
Mesh cullMesh(const Mesh& input_mesh, const Camera& camera);

}  // namespace sim