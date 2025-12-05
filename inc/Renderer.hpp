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

#include <iostream>

#include "Vector.hpp"
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
 * @brief Internal type for pixel position (framebuffer indexing only).
 */
typedef struct {
  int_t x;
  int_t y;
} PixelPosition;

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
  int_t m_width, m_height;
  std::vector<Color> m_pixel_colors;

public:
  FrameBuffer(int_t width, int_t height);

  Color& operator[](PixelPosition px);
  const Color& operator[](PixelPosition px) const;

  void clear(const Color& color = black);

  int_t getWidth() const;
  int_t getHeight() const;
  const void* getRawData() const;
};

#pragma endregion

#pragma region Window

class Window
{
private:
  SDL_Window* m_window;
  SDL_Renderer* m_SDLRenderer;
  SDL_Texture* m_Texture;
  int_t m_width, m_height;

public:
  Window(int_t w, int_t h);

  ~Window();

  void display(const FrameBuffer& frameBuffer);
  bool shouldClose();
};

#pragma endregion

#pragma region Renderer

class Renderer
{
private:
  // FrameBuffer and SDL members (formerly Window)
  FrameBuffer m_buffer;
  SDL_Window* m_window;
  SDL_Renderer* m_SDLRenderer;
  SDL_Texture* m_texture;

  // Helper members for coordinate system and scaling
  real_t m_scale;           // Pixels per unit (e.g., 100 pixels per 1 unit)
  PixelPosition m_center;   // Center of the viewport in pixels
  int_t m_grid_spacing;     // d: grid spacing in pixels (1/scale)
  int_t m_width, m_height;  // window width & height

  /**
   * @brief Convert world coordinates (Vector) to pixel coordinates.
   * @param worldPos 2D world position vector.
   * @return Pixel position for framebuffer access.
   * @throws std::invalid_argument if vector is not 2D.
   */
  PixelPosition worldToPixel(const Vector& worldPos) const;

public:
  /**
   * @brief Construct a Renderer with given width, height and scaling.
   *
   * @param width Window width in pixels.
   * @param height Window height in pixels.
   * @param scale Pixels per unit (default: 100).
   */
  Renderer(int_t width, int_t height, real_t scale = 100.0);

  /**
   * @brief Destructor.
   */
  ~Renderer();

  void drawLine(const Vector& start, const Vector& end, uint_t widthPixels, const Color& color);
  void drawCircle(const Vector& center, real_t radius, const Color& color);
  void drawPolygon(const std::vector<Vector>& vertices, uint_t widthPixels, const Color& color);
  void drawArrow(const Vector& start, const Vector& end, uint_t widthPixels, const Color& color);
  void drawGrid(uint_t widthPixels = 1, const Color& color = gray);
  void drawWorldFrame(uint_t widthPixels = 2);

  /**
   * @brief Display the rendered frame and handle events until window closes.
   *
   * @details This is a blocking call that shows the window and runs the
   * event loop, processing SDL events and checking if the window should close.
   * The frame rate is capped at ~60 FPS. Call with no arguments.
   */
  void display();

  void clear(const Color& color = black) { m_buffer.clear(color); }
};

#pragma endregion

}  // namespace sim