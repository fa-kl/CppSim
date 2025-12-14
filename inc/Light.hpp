/*****************************************************************************************
 * @file: Light.hpp
 *
 * @brief: Light sources for 3D rendering.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include "Color.hpp"
#include "Vector.hpp"
#include "types.hpp"

namespace sim
{

/**
 * @brief Types of light sources.
 */
enum class LightType
{
    AMBIENT,      ///< Ambient light (no direction)
    DIRECTIONAL,  ///< Directional light (like sun)
    POINT         ///< Point light (has position)
};

/**
 * @brief A light source for illuminating 3D scenes.
 */
class Light
{
  protected:
    LightType m_type;
    Color m_color;
    real_t m_intensity;
    Vector m_direction;  // For directional lights
    Vector m_position;   // For point lights

  public:
    /**
     * @brief Create an ambient light.
     * @param color Light color.
     * @param intensity Light intensity (0-1).
     */
    static Light ambient(const Color& color, real_t intensity = 0.2);

    /**
     * @brief Create a directional light.
     * @param direction Light direction (will be normalized).
     * @param color Light color.
     * @param intensity Light intensity (0-1).
     */
    static Light directional(const Vector& direction, const Color& color, real_t intensity = 1.0);

    /**
     * @brief Create a point light.
     * @param position Light position in world space.
     * @param color Light color.
     * @param intensity Light intensity (0-1).
     */
    static Light point(const Vector& position, const Color& color, real_t intensity = 1.0);

    /**
     * @brief Get the light type.
     */
    LightType getType() const;

    /**
     * @brief Get the light color.
     */
    Color getColor() const;

    /**
     * @brief Get the light intensity.
     */
    real_t getIntensity() const;

    /**
     * @brief Get the light direction (for directional lights).
     */
    Vector getDirection() const;

    /**
     * @brief Get the light position (for point lights).
     */
    Vector getPosition() const;

  private:
    Light(LightType type);
};

}  // namespace sim
