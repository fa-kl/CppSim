/*****************************************************************************************
 * @file: Vertex.cpp
 *
 * @brief: Vertex struct implementations for rendering.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Vertex.hpp"

#include <stdexcept>

namespace sim
{

Vertex::Vertex(const Vector& p, const Color& c, const Vector& n) : position(p), color(c), normal(normalize(n))
{
    if (position.length() != 3)
    {
        throw std::invalid_argument("Position must be 3-dimensional");
    }
    if (normal.length() != 3)
    {
        throw std::invalid_argument("Normal must be 3-dimensional");
    }
}

}  // namespace sim