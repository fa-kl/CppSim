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

#include <algorithm>

#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

struct Color
{
    uint8_t red;
    uint8_t green;
    uint8_t blue;
    uint8_t opacity;

    /* Constructors */
    Color();
    Color(uint8_t r, uint8_t g, uint8_t b, uint8_t o);
    Color(real_t r, real_t g, real_t b, real_t o);

    /* Arithmetic operators */
    Color operator+(const Color& other) const;
    Color operator-(const Color& other) const;
    Color operator*(real_t scalar) const;
    Color operator/(real_t scalar) const;
    Color& operator+=(const Color& other);
    Color& operator-=(const Color& other);
    Color& operator*=(real_t scalar);
    Color& operator/=(real_t scalar);

    /* Static color constants */
    static Color White();
    static Color Black();
    static Color Transparent();
    static Color Red();
    static Color Green();
    static Color Blue();
    static Color Yellow();
    static Color Cyan();
    static Color Magenta();
    static Color Orange();
    static Color Purple();
    static Color Pink();
    static Color Gray();
    static Color Silver();
    static Color LightRed();
    static Color LightGreen();
    static Color LightBlue();
    static Color LightYellow();
    static Color LightCyan();
    static Color LightMagenta();
    static Color LightOrange();
    static Color LightPurple();
    static Color LightPink();
    static Color LightGray();
    static Color DarkRed();
    static Color DarkGreen();
    static Color DarkBlue();
    static Color DarkYellow();
    static Color DarkCyan();
    static Color DarkMagenta();
    static Color DarkOrange();
    static Color DarkPurple();
    static Color DarkPink();
    static Color DarkGray();
};

/* Non-member operators*/
Color operator*(real_t scalar, const Color& color);

/* Non-member functions */
Color mean(std::vector<Color> colors);

/**
 * @brief Linear interpolation between two colors.
 * @param c1 First color.
 * @param c2 Second color.
 * @param t Interpolation parameter (0.0 = c1, 1.0 = c2).
 * @return Interpolated color.
 */
Color interpolate(const Color& c1, const Color& c2, real_t t);

/**
 * @brief Barycentric interpolation between three colors.
 * @param c1 First color.
 * @param c2 Second color.
 * @param c3 Third color.
 * @param weights Barycentric weights (3D vector: w, u, v).
 * @return Interpolated color.
 */
Color interpolate(const Color& c1, const Color& c2, const Color& c3, const Vector& weights);

}  // namespace sim
