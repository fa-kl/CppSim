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
#include "FrameBuffer.hpp"
#include "Matrix.hpp"
#include "Mesh.hpp"
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
    void drawLine(PixelCoordinate px1, PixelCoordinate px2, const Color& color);

    /**
     * @brief Fill a triangle using scanline rasterization.
     */
    void fillTriangle(PixelCoordinate px1, PixelCoordinate px2, PixelCoordinate px3, const Color& color);

    /**
     * @brief Convert NDC coordinates to screen coordinates.
     */
    PixelCoordinate ndcToScreen(const Vector& ndc) const;

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

}  // namespace sim