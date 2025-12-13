/*****************************************************************************************
 * @file: Color.hpp
 *
 * @brief: Color type for rendering.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include "types.hpp"

namespace sim
{

struct Color
{
    uint8_t red;
    uint8_t green;
    uint8_t blue;
    uint8_t opacity;

    Color();

    Color(uint8_t r, uint8_t g, uint8_t b, uint8_t o);

    Color(real_t r, real_t g, real_t b, real_t o);
};

Color mean(std::vector<Color> colors);

extern const Color WHITE;
extern const Color BLACK;
extern const Color TRANSPARENT;
extern const Color RED;
extern const Color GREEN;
extern const Color BLUE;
extern const Color YELLOW;
extern const Color CYAN;
extern const Color MAGENTA;
extern const Color ORANGE;
extern const Color PURPLE;
extern const Color PINK;
extern const Color GRAY;
extern const Color SILVER;
extern const Color LIGHT_RED;
extern const Color LIGHT_GREEN;
extern const Color LIGHT_BLUE;
extern const Color LIGHT_YELLOW;
extern const Color LIGHT_CYAN;
extern const Color LIGHT_MAGENTA;
extern const Color LIGHT_ORANGE;
extern const Color LIGHT_PURPLE;
extern const Color LIGHT_PINK;
extern const Color LIGHT_GRAY;
extern const Color DARK_RED;
extern const Color DARK_GREEN;
extern const Color DARK_BLUE;
extern const Color DARK_YELLOW;
extern const Color DARK_CYAN;
extern const Color DARK_MAGENTA;
extern const Color DARK_ORANGE;
extern const Color DARK_PURPLE;
extern const Color DARK_PINK;
extern const Color DARK_GRAY;

}  // namespace sim
