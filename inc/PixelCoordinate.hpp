/*****************************************************************************************
 * @file: PixelCoordinate.hpp
 *
 * @brief: This file provides a container for pixel coordinates.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include <cmath>
#include <stdexcept>

#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

struct PixelCoordinate
{
    int_t x;
    int_t y;

    /* Constructors */

    inline PixelCoordinate() : x(0), y(0) {}

    inline PixelCoordinate(int_t px, int_t py) : x(px), y(py) {}

    inline PixelCoordinate(const Vector& vec)
    {
        if (vec.length() != 2)
        {
            throw std::invalid_argument("Vector must be 2-dimensional to represent a pixel coordinate");
        }
        x = static_cast<int_t>(std::round(vec[0]));
        y = static_cast<int_t>(std::round(vec[1]));
    }
};

}  // namespace sim