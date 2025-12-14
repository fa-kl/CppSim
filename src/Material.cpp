/*****************************************************************************************
 * @file: Material.cpp
 *
 * @brief: Implementation of Material struct.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Material.hpp"

#include <stdexcept>

namespace sim
{

Material::Material(real_t rho, Color c) : density{rho}, color{c}
{
    if (density <= 0)
    {
        throw std::invalid_argument("Density must be greater than zero");
    }
}

}  // namespace sim
