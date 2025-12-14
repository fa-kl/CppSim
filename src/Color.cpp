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

Color interpolate(const Color& c1, const Color& c2, real_t t)
{
    // Clamp t to [0, 1]
    t = std::clamp(t, 0.0, 1.0);
    real_t s = 1.0 - t;

    real_t r = s * static_cast<real_t>(c1.red) + t * static_cast<real_t>(c2.red);
    real_t g = s * static_cast<real_t>(c1.green) + t * static_cast<real_t>(c2.green);
    real_t b = s * static_cast<real_t>(c1.blue) + t * static_cast<real_t>(c2.blue);
    real_t a = s * static_cast<real_t>(c1.opacity) + t * static_cast<real_t>(c2.opacity);

    Color result;
    result.red = static_cast<uint8_t>(std::clamp(r, 0.0, 255.0));
    result.green = static_cast<uint8_t>(std::clamp(g, 0.0, 255.0));
    result.blue = static_cast<uint8_t>(std::clamp(b, 0.0, 255.0));
    result.opacity = static_cast<uint8_t>(std::clamp(a, 0.0, 255.0));

    return result;
}

Color interpolate(const Color& c1, const Color& c2, const Color& c3, const Vector& weights)
{
    if (weights.length() != 3)
    {
        throw std::invalid_argument("Weights vector must have exactly 3 elements for barycentric interpolation");
    }

    real_t w = weights[0];
    real_t u = weights[1];
    real_t v = weights[2];

    real_t r = w * static_cast<real_t>(c1.red) + u * static_cast<real_t>(c2.red) + v * static_cast<real_t>(c3.red);
    real_t g =
        w * static_cast<real_t>(c1.green) + u * static_cast<real_t>(c2.green) + v * static_cast<real_t>(c3.green);
    real_t b = w * static_cast<real_t>(c1.blue) + u * static_cast<real_t>(c2.blue) + v * static_cast<real_t>(c3.blue);
    real_t a =
        w * static_cast<real_t>(c1.opacity) + u * static_cast<real_t>(c2.opacity) + v * static_cast<real_t>(c3.opacity);

    Color result;
    result.red = static_cast<uint8_t>(std::clamp(r, 0.0, 255.0));
    result.green = static_cast<uint8_t>(std::clamp(g, 0.0, 255.0));
    result.blue = static_cast<uint8_t>(std::clamp(b, 0.0, 255.0));
    result.opacity = static_cast<uint8_t>(std::clamp(a, 0.0, 255.0));

    return result;
}

Color Color::operator+(const Color& other) const
{
    return Color{static_cast<uint8_t>(std::min(255, static_cast<int>(red) + static_cast<int>(other.red))),
                 static_cast<uint8_t>(std::min(255, static_cast<int>(green) + static_cast<int>(other.green))),
                 static_cast<uint8_t>(std::min(255, static_cast<int>(blue) + static_cast<int>(other.blue))),
                 static_cast<uint8_t>(std::min(255, static_cast<int>(opacity) + static_cast<int>(other.opacity)))};
}

Color Color::operator-(const Color& other) const
{
    return Color{static_cast<uint8_t>(std::max(0, static_cast<int>(red) - static_cast<int>(other.red))),
                 static_cast<uint8_t>(std::max(0, static_cast<int>(green) - static_cast<int>(other.green))),
                 static_cast<uint8_t>(std::max(0, static_cast<int>(blue) - static_cast<int>(other.blue))),
                 static_cast<uint8_t>(std::max(0, static_cast<int>(opacity) - static_cast<int>(other.opacity)))};
}

Color Color::operator*(real_t scalar) const
{
    return Color{static_cast<uint8_t>(std::clamp(static_cast<int>(red * scalar), 0, 255)),
                 static_cast<uint8_t>(std::clamp(static_cast<int>(green * scalar), 0, 255)),
                 static_cast<uint8_t>(std::clamp(static_cast<int>(blue * scalar), 0, 255)),
                 opacity};
}

Color Color::operator/(real_t scalar) const
{
    if (scalar == 0.0)
    {
        return Color::Black();
    }
    return *this * (1.0 / scalar);
}

Color& Color::operator+=(const Color& other)
{
    *this = *this + other;
    return *this;
}

Color& Color::operator-=(const Color& other)
{
    *this = *this - other;
    return *this;
}

Color& Color::operator*=(real_t scalar)
{
    *this = *this * scalar;
    return *this;
}

Color& Color::operator/=(real_t scalar)
{
    *this = *this / scalar;
    return *this;
}

Color operator*(real_t scalar, const Color& color)
{
    return color * scalar;
}

Color Color::White()
{
    return Color((uint8_t)255, (uint8_t)255, (uint8_t)255, (uint8_t)255);
}

Color Color::Black()
{
    return Color((uint8_t)0, (uint8_t)0, (uint8_t)0, (uint8_t)255);
}

Color Color::Transparent()
{
    return Color((uint8_t)0, (uint8_t)0, (uint8_t)0, (uint8_t)0);
}

Color Color::Red()
{
    return Color((uint8_t)255, (uint8_t)0, (uint8_t)0, (uint8_t)255);
}

Color Color::Green()
{
    return Color((uint8_t)0, (uint8_t)255, (uint8_t)0, (uint8_t)255);
}

Color Color::Blue()
{
    return Color((uint8_t)0, (uint8_t)0, (uint8_t)255, (uint8_t)255);
}

Color Color::Yellow()
{
    return Color((uint8_t)255, (uint8_t)255, (uint8_t)0, (uint8_t)255);
}

Color Color::Cyan()
{
    return Color((uint8_t)0, (uint8_t)255, (uint8_t)255, (uint8_t)255);
}

Color Color::Magenta()
{
    return Color((uint8_t)255, (uint8_t)0, (uint8_t)255, (uint8_t)255);
}

Color Color::Orange()
{
    return Color((uint8_t)255, (uint8_t)165, (uint8_t)0, (uint8_t)255);
}

Color Color::Purple()
{
    return Color((uint8_t)128, (uint8_t)0, (uint8_t)128, (uint8_t)255);
}

Color Color::Pink()
{
    return Color((uint8_t)255, (uint8_t)192, (uint8_t)203, (uint8_t)255);
}

Color Color::Gray()
{
    return Color((uint8_t)128, (uint8_t)128, (uint8_t)128, (uint8_t)255);
}

Color Color::Silver()
{
    return Color((uint8_t)192, (uint8_t)192, (uint8_t)192, (uint8_t)255);
}

Color Color::LightRed()
{
    return Color((uint8_t)255, (uint8_t)128, (uint8_t)128, (uint8_t)255);
}

Color Color::LightGreen()
{
    return Color((uint8_t)144, (uint8_t)238, (uint8_t)144, (uint8_t)255);
}

Color Color::LightBlue()
{
    return Color((uint8_t)173, (uint8_t)216, (uint8_t)230, (uint8_t)255);
}

Color Color::LightYellow()
{
    return Color((uint8_t)255, (uint8_t)255, (uint8_t)224, (uint8_t)255);
}

Color Color::LightCyan()
{
    return Color((uint8_t)224, (uint8_t)255, (uint8_t)255, (uint8_t)255);
}

Color Color::LightMagenta()
{
    return Color((uint8_t)255, (uint8_t)128, (uint8_t)255, (uint8_t)255);
}

Color Color::LightOrange()
{
    return Color((uint8_t)255, (uint8_t)200, (uint8_t)124, (uint8_t)255);
}

Color Color::LightPurple()
{
    return Color((uint8_t)200, (uint8_t)162, (uint8_t)200, (uint8_t)255);
}

Color Color::LightPink()
{
    return Color((uint8_t)255, (uint8_t)220, (uint8_t)225, (uint8_t)255);
}

Color Color::LightGray()
{
    return Color((uint8_t)211, (uint8_t)211, (uint8_t)211, (uint8_t)255);
}

Color Color::DarkRed()
{
    return Color((uint8_t)139, (uint8_t)0, (uint8_t)0, (uint8_t)255);
}

Color Color::DarkGreen()
{
    return Color((uint8_t)0, (uint8_t)100, (uint8_t)0, (uint8_t)255);
}

Color Color::DarkBlue()
{
    return Color((uint8_t)0, (uint8_t)0, (uint8_t)139, (uint8_t)255);
}

Color Color::DarkYellow()
{
    return Color((uint8_t)204, (uint8_t)204, (uint8_t)0, (uint8_t)255);
}

Color Color::DarkCyan()
{
    return Color((uint8_t)0, (uint8_t)139, (uint8_t)139, (uint8_t)255);
}

Color Color::DarkMagenta()
{
    return Color((uint8_t)139, (uint8_t)0, (uint8_t)139, (uint8_t)255);
}

Color Color::DarkOrange()
{
    return Color((uint8_t)255, (uint8_t)140, (uint8_t)0, (uint8_t)255);
}

Color Color::DarkPurple()
{
    return Color((uint8_t)75, (uint8_t)0, (uint8_t)130, (uint8_t)255);
}

Color Color::DarkPink()
{
    return Color((uint8_t)231, (uint8_t)84, (uint8_t)128, (uint8_t)255);
}

Color Color::DarkGray()
{
    return Color((uint8_t)169, (uint8_t)169, (uint8_t)169, (uint8_t)255);
}

}  // namespace sim