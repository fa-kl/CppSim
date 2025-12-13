/*****************************************************************************************
 * @file: FrameBuffer.cpp
 *
 * @brief: This file implements a 2D frame buffer class.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "FrameBuffer.hpp"

#include <stdexcept>

namespace sim
{

bool FrameBuffer::isInFrame(const PixelCoordinate& px) const
{
    return 0 <= px.x && px.x < m_width && 0 <= px.y && px.y < m_height;
}

FrameBuffer::FrameBuffer(int_t width, int_t height)
    : m_width(width), m_height(height), m_pixel_colors(static_cast<uint_t>(width * height))
{
    if (width <= 0 || height <= 0)
    {
        throw std::invalid_argument("Width and height must both be greater than zero");
    }
}

Color& FrameBuffer::operator[](const PixelCoordinate& px)
{
    if (isInFrame(px))
    {
        return m_pixel_colors[static_cast<uint_t>(px.y * m_width + px.x)];
    }
    std::string pixel = "(" + std::to_string(px.x) + ", " + std::to_string(px.y) + ")";
    std::string window_size = "(" + std::to_string(m_width) + ", " + std::to_string(m_height) + ")";
    std::string msg = "Pixel " + pixel + "out of range for window of size " + window_size + ".";
    throw std::out_of_range(msg);
}

const Color& FrameBuffer::operator[](const PixelCoordinate& px) const
{
    if (isInFrame(px))
    {
        return m_pixel_colors[static_cast<uint_t>(px.y * m_width + px.x)];
    }
    std::string pixel = "(" + std::to_string(px.x) + ", " + std::to_string(px.y) + ")";
    std::string window_size = "(" + std::to_string(m_width) + ", " + std::to_string(m_height) + ")";
    std::string msg = "Pixel " + pixel + "out of range for window of size " + window_size + ".";
    throw std::out_of_range(msg);
}

void FrameBuffer::clear(const Color& color)
{
    std::fill(m_pixel_colors.begin(), m_pixel_colors.end(), color);
}

int_t FrameBuffer::getWidth() const
{
    return m_width;
}

int_t FrameBuffer::getHeight() const
{
    return m_height;
}

const void* FrameBuffer::getRawData() const
{
    return m_pixel_colors.data();
}

}  // namespace sim