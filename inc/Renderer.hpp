/*****************************************************************************************
 * @file: Renderer.hpp
 *
 * @brief: This file provides all types, classes and functions for rendering.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#pragma once

#include "types.hpp"

#ifdef __linux__
#include <SDL2/SDL.h>
#else
#include <SDL.h>
#endif

namespace sim
{

#pragma region PixelPosition

/**
 * @brief A type for a pixel position.
 */
typedef struct {
  uint_t x;
  uint_t y;
} PixelPosition;

/**
 * @brief Add two pixel positions.
 */
PixelPosition operator+(PixelPosition px1, const PixelPosition& px2);

/**
 * @brief Subtract two pixel positions.
 */
PixelPosition operator-(PixelPosition px1, const PixelPosition& px2);

/**
 * @brief Multiply a pixel position with a scalar integer.
 */
PixelPosition operator*(PixelPosition px, const int_t& k);

/**
 * @brief Multiply a pixel position with a scalar integer.
 */
PixelPosition operator*(const int_t& k, PixelPosition px);

#pragma endregion

#pragma region Color

/**
 * @brief A type for storing color values in RGBA format.
 */
typedef struct {
  uint8_t red;
  uint8_t green;
  uint8_t blue;
  uint8_t opacity;
} Color;

/**
 * @brief Create a color using real values from the interval [0; 1].
 */
Color color(real_t red, real_t green, real_t blue, real_t opacity);

extern const Color white;
extern const Color black;
extern const Color transparent;
extern const Color red;
extern const Color green;
extern const Color blue;
extern const Color yellow;
extern const Color cyan;
extern const Color magenta;
extern const Color orange;
extern const Color purple;
extern const Color pink;
extern const Color lightRed;
extern const Color lightGreen;
extern const Color lightBlue;
extern const Color lightYellow;
extern const Color lightCyan;
extern const Color lightMagenta;
extern const Color lightOrange;
extern const Color lightPurple;
extern const Color lightPink;
extern const Color lightGray;
extern const Color darkRed;
extern const Color darkGreen;
extern const Color darkBlue;
extern const Color darkYellow;
extern const Color darkCyan;
extern const Color darkMagenta;
extern const Color darkOrange;
extern const Color darkPurple;
extern const Color darkPink;
extern const Color darkGray;
extern const Color gray;
extern const Color silver;
extern const Color dimGray;
extern const Color redAlpha;
extern const Color greenAlpha;
extern const Color blueAlpha;
extern const Color yellowAlpha;
extern const Color cyanAlpha;
extern const Color magentaAlpha;
extern const Color whiteAlpha;
extern const Color blackAlpha;

#pragma endregion

#pragma region FrameBuffer

class FrameBuffer
{
protected:
  uint_t m_width, m_height;
  std::vector<Color> m_pixel_colors;

public:
  FrameBuffer(uint_t width, uint_t height);

  Color& operator[](PixelPosition px);
  const Color& operator[](PixelPosition px) const;

  void clear(const Color& color = black);

  uint_t getWidth() const;
  uint_t getHeight() const;
  const void* getRawData() const;
};

#pragma endregion

#pragma region Renderer

class Renderer
{
protected:
  FrameBuffer& m_buffer;

public:
  Renderer(FrameBuffer& fb);

  void drawLine(const PixelPosition& start, const PixelPosition& end, uint_t width, const Color& color);
  void drawCircle(const PixelPosition& center, const real_t& radius, const Color& color);
  void drawPolygon(const std::vector<PixelPosition>& vertices, const uint_t& lineWidth, const Color& color);
};

#pragma endregion

#pragma region Window

class Window
{
private:
  SDL_Window* window;
  SDL_Renderer* sdlRenderer;
  SDL_Texture* texture;
  uint_t width, height;

public:
  Window(uint_t w, uint_t h);

  ~Window();

  void display(const FrameBuffer& frameBuffer);
  bool shouldClose();
};

#pragma endregion

}  // namespace sim