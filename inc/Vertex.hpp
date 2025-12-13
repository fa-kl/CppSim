/*****************************************************************************************
 * @file: Vertex.hpp
 *
 * @brief: Vertex type for rendering.
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

struct Vertex {
  Vector point;
  Color color;

  Vertex(Vector p, Color c);
};

}  // namespace sim