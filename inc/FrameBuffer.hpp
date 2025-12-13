/*****************************************************************************************
 * @file: FrameBuffer.hpp
 *
 * @brief: This file provides a 2D frame buffer class.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include <vector>

#include "Color.hpp"
#include "types.hpp"

namespace sim
{

struct PixelCoordinate
{
    int_t x;
    int_t y;
};

class FrameBuffer
{
  protected:
    int_t m_width;
    int_t m_height;
    std::vector<Color> m_pixel_colors;

    bool isInFrame(const PixelCoordinate& px) const;

  public:
    FrameBuffer(int_t width, int_t height);

    Color& operator[](const PixelCoordinate& px);
    const Color& operator[](const PixelCoordinate& px) const;

    void clear(const Color& color = BLACK);

    int_t getWidth() const;
    int_t getHeight() const;
    const void* getRawData() const;
};

}  // namespace sim