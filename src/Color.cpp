/*****************************************************************************************
 * @file: Color.cpp
 *
 * @brief: Color type implementations for rendering.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Color.hpp"

#include <cmath>
#include <stdexcept>

namespace sim
{

Color::Color()
{
    red = 0;
    blue = 0;
    green = 0;
    opacity = 0;
}

Color::Color(uint8_t r, uint8_t g, uint8_t b, uint8_t o)
{
    red = r;
    green = g;
    blue = b;
    opacity = o;
}

Color::Color(real_t r, real_t g, real_t b, real_t o)
{
    if (r < 0 || r > 1)
    {
        throw std::invalid_argument("All values must be within [0; 1].");
    }
    if (g < 0 || g > 1)
    {
        throw std::invalid_argument("All values must be within [0; 1].");
    }
    if (b < 0 || b > 1)
    {
        throw std::invalid_argument("All values must be within [0; 1].");
    }
    if (o < 0 || o > 1)
    {
        throw std::invalid_argument("All values must be within [0; 1].");
    }
    red = static_cast<uint8_t>(round(255.0 * r));
    green = static_cast<uint8_t>(round(255.0 * g));
    blue = static_cast<uint8_t>(round(255.0 * b));
    opacity = static_cast<uint8_t>(round(255.0 * o));
}

Color mean(std::vector<Color> colors)
{
    real_t mean_red = 0.0;
    real_t mean_green = 0.0;
    real_t mean_blue = 0.0;
    real_t mean_opacity = 0.0;
    for (const Color& col : colors)
    {
        mean_red += static_cast<real_t>(col.red) / 255.0;
        mean_green += static_cast<real_t>(col.green) / 255.0;
        mean_blue += static_cast<real_t>(col.blue) / 255.0;
        mean_opacity += static_cast<real_t>(col.opacity) / 255.0;
    }
    real_t n = static_cast<real_t>(colors.size());
    return {mean_red / n, mean_green / n, mean_blue / n, mean_opacity / n};
}

const Color WHITE = {(uint8_t)255, (uint8_t)255, (uint8_t)255, (uint8_t)255};
const Color BLACK = {(uint8_t)0, (uint8_t)0, (uint8_t)0, (uint8_t)255};
const Color TRANSPARENT = {(uint8_t)0, (uint8_t)0, (uint8_t)0, (uint8_t)0};
const Color RED = {(uint8_t)255, (uint8_t)0, (uint8_t)0, (uint8_t)255};
const Color GREEN = {(uint8_t)0, (uint8_t)255, (uint8_t)0, (uint8_t)255};
const Color BLUE = {(uint8_t)0, (uint8_t)0, (uint8_t)255, (uint8_t)255};
const Color YELLOW = {(uint8_t)255, (uint8_t)255, (uint8_t)0, (uint8_t)255};
const Color CYAN = {(uint8_t)0, (uint8_t)255, (uint8_t)255, (uint8_t)255};
const Color MAGENTA = {(uint8_t)255, (uint8_t)0, (uint8_t)255, (uint8_t)255};
const Color ORANGE = {(uint8_t)255, (uint8_t)165, (uint8_t)0, (uint8_t)255};
const Color PURPLE = {(uint8_t)128, (uint8_t)0, (uint8_t)128, (uint8_t)255};
const Color PINK = {(uint8_t)255, (uint8_t)192, (uint8_t)203, (uint8_t)255};
const Color GRAY = {(uint8_t)128, (uint8_t)128, (uint8_t)128, (uint8_t)255};
const Color SILVER = {(uint8_t)192, (uint8_t)192, (uint8_t)192, (uint8_t)255};
const Color LIGHT_RED = {(uint8_t)255, (uint8_t)128, (uint8_t)128, (uint8_t)255};
const Color LIGHT_GREEN = {(uint8_t)144, (uint8_t)238, (uint8_t)144, (uint8_t)255};
const Color LIGHT_BLUE = {(uint8_t)173, (uint8_t)216, (uint8_t)230, (uint8_t)255};
const Color LIGHT_YELLOW = {(uint8_t)255, (uint8_t)255, (uint8_t)224, (uint8_t)255};
const Color LIGHT_CYAN = {(uint8_t)224, (uint8_t)255, (uint8_t)255, (uint8_t)255};
const Color LIGHT_MAGENTA = {(uint8_t)255, (uint8_t)128, (uint8_t)255, (uint8_t)255};
const Color LIGHT_ORANGE = {(uint8_t)255, (uint8_t)200, (uint8_t)124, (uint8_t)255};
const Color LIGHT_PURPLE = {(uint8_t)200, (uint8_t)162, (uint8_t)200, (uint8_t)255};
const Color LIGHT_PINK = {(uint8_t)255, (uint8_t)220, (uint8_t)225, (uint8_t)255};
const Color LIGHT_GRAY = {(uint8_t)211, (uint8_t)211, (uint8_t)211, (uint8_t)255};
const Color DARK_RED = {(uint8_t)139, (uint8_t)0, (uint8_t)0, (uint8_t)255};
const Color DARK_GREEN = {(uint8_t)0, (uint8_t)100, (uint8_t)0, (uint8_t)255};
const Color DARK_BLUE = {(uint8_t)0, (uint8_t)0, (uint8_t)139, (uint8_t)255};
const Color DARK_YELLOW = {(uint8_t)204, (uint8_t)204, (uint8_t)0, (uint8_t)255};
const Color DARK_CYAN = {(uint8_t)0, (uint8_t)139, (uint8_t)139, (uint8_t)255};
const Color DARK_MAGENTA = {(uint8_t)139, (uint8_t)0, (uint8_t)139, (uint8_t)255};
const Color DARK_ORANGE = {(uint8_t)255, (uint8_t)140, (uint8_t)0, (uint8_t)255};
const Color DARK_PURPLE = {(uint8_t)75, (uint8_t)0, (uint8_t)130, (uint8_t)255};
const Color DARK_PINK = {(uint8_t)231, (uint8_t)84, (uint8_t)128, (uint8_t)255};
const Color DARK_GRAY = {(uint8_t)169, (uint8_t)169, (uint8_t)169, (uint8_t)255};

}  // namespace sim