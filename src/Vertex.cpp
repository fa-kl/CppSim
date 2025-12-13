/*****************************************************************************************
 * @file: Vertex.cpp
 *
 * @brief: Vertex type implementations for rendering.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Vertex.hpp"

#include <stdexcept>

namespace sim
{

Vertex::Vertex(Vector p, Color c) : point(p), color(c)
{
  if (point.length() != 3) {
    throw std::invalid_argument("Point must be 3-dimensional");
  }
}

}  // namespace sim