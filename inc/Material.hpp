/*****************************************************************************************
 * @file: Material.hpp
 *
 * @brief: Material properties for lighting calculations and physics.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include "Color.hpp"
#include "types.hpp"

namespace sim
{

/**
 * @brief Material properties combining physics and rendering.
 */
struct Material
{
    real_t density;
    Color color;

    Material(real_t density, Color color);
};

}  // namespace sim