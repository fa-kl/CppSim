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

#include <array>
#include <iostream>
#include <vector>

#include "Camera.hpp"
#include "Color.hpp"
#include "Matrix.hpp"
#include "Rotation.hpp"
#include "Shapes.hpp"
#include "Transform.hpp"
#include "Vector.hpp"
#include "Vertex.hpp"
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

#pragma region FrameBuffer

class FrameBuffer
{
protected:
  int_t m_width, m_height;
  std::vector<Color> m_pixel_colors;

  bool isInFrame(const PixelPosition& px) const;

public:
  FrameBuffer(int_t width, int_t height);

  Color& operator[](PixelPosition px);
  const Color& operator[](PixelPosition px) const;

  void clear(const Color& color = BLACK);

  int_t getWidth() const;
  int_t getHeight() const;
  const void* getRawData() const;
};

#pragma endregion

#pragma region MeshTransform

/**
 * @brief Apply a transformation to a mesh.
 */
Mesh operator*(const Transform<3>& transform, const Mesh& mesh);

/**
 * @brief Cull triangles that are completely outside the camera's view frustum.
 */
Mesh cullMesh(const Mesh& input_mesh, const Camera& camera);

#pragma endregion

#pragma region Renderer

/**
 * @brief 3D Renderer using SDL2 and a framebuffer.
 */
class Renderer
{
protected:
  SDL_Window* m_window;
  SDL_Renderer* m_sdl_renderer;
  SDL_Texture* m_texture;
  FrameBuffer m_framebuffer;
  int_t m_width, m_height;
  bool m_should_close;

  /**
   * @brief Draw a line using Bresenham's algorithm.
   */
  void drawLine(int x0, int y0, int x1, int y1, const Color& color);

  /**
   * @brief Fill a triangle using scanline rasterization.
   */
  void fillTriangle(int x0, int y0, int x1, int y1, int x2, int y2, const Color& color);

  /**
   * @brief Convert NDC coordinates to screen coordinates.
   */
  std::pair<int, int> ndcToScreen(const Vector& ndc) const;

public:
  /**
   * @brief Construct a renderer with specified window size.
   */
  Renderer(int_t width, int_t height, const std::string& title = "3D Renderer");

  /**
   * @brief Destructor - cleanup SDL resources.
   */
  ~Renderer();

  /**
   * @brief Clear the framebuffer with a color.
   */
  void clear(const Color& color = BLACK);

  /**
   * @brief Render a mesh to the framebuffer.
   */
  void renderMesh(const Mesh& mesh, const Camera& camera, const Color& line_color = WHITE);

  /**
   * @brief Draw coordinate axes at origin.
   */
  void drawAxes(const Camera& camera, real_t length = 1.0);

  /**
   * @brief Present the framebuffer to the screen.
   */
  void present();

  /**
   * @brief Get the framebuffer for direct manipulation.
   */
  FrameBuffer& getFrameBuffer();

  /**
   * @brief Get window width.
   */
  int_t getWidth() const;

  /**
   * @brief Get window height.
   */
  int_t getHeight() const;

  /**
   * @brief Check if window should close.
   */
  bool shouldClose();

  /**
   * @brief Poll events and update window state.
   */
  void pollEvents();
};

#pragma endregion

}  // namespace sim