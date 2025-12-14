/*****************************************************************************************
 * @file: Light.cpp
 *
 * @brief: Implementation of Light class.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Light.hpp"

#include <cmath>

namespace sim
{

Light::Light(LightType type)
    : m_type{type}, m_color{Color::White()}, m_intensity{1.0}, m_direction{0, 0, -1}, m_position{0, 0, 0}
{
}

Light Light::ambient(const Color& color, real_t intensity)
{
    Light light(LightType::AMBIENT);
    light.m_color = color;
    light.m_intensity = intensity;
    return light;
}

Light Light::directional(const Vector& direction, const Color& color, real_t intensity)
{
    Light light(LightType::DIRECTIONAL);
    light.m_color = color;
    light.m_intensity = intensity;

    // Normalize direction
    real_t length = std::sqrt(direction[1] * direction[1] + direction[2] * direction[2] + direction[3] * direction[3]);
    if (length > 0)
    {
        light.m_direction = Vector{direction[1] / length, direction[2] / length, direction[3] / length};
    }

    return light;
}

Light Light::point(const Vector& position, const Color& color, real_t intensity)
{
    Light light(LightType::POINT);
    light.m_color = color;
    light.m_intensity = intensity;
    light.m_position = position;
    return light;
}

LightType Light::getType() const
{
    return m_type;
}

Color Light::getColor() const
{
    return m_color;
}

real_t Light::getIntensity() const
{
    return m_intensity;
}

Vector Light::getDirection() const
{
    return m_direction;
}

Vector Light::getPosition() const
{
    return m_position;
}

}  // namespace sim
