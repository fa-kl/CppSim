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
#include "RigidBody.hpp"
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
    Camera m_camera;
    bool m_should_close;

    /**
     * @brief Draw a 1 pixel wide line between two world coordinates.
     */
    void drawLine(const Vector& v1, const Vector& v2, const Color& color);

    /**
     * @brief Fill a triangle based on it's world vertices.
     *
     * @details This method not only projects the vertices given in world coordinates
     * onto the screen but also applies the correct colors to each pixel based on the
     * vertex colors.
     */
    void fillTriangle(const Vertex& v1, const Vertex& v2, const Vertex& v3);

  public:
    /**
     * @brief Construct a renderer with specified window size and camera.
     * @param width Window width in pixels.
     * @param height Window height in pixels.
     * @param camera Camera for rendering (will be fixed after construction).
     * @param title Window title.
     */
    Renderer(const Camera& camera, const std::string& title = "3D Renderer");

    /**
     * @brief Destructor - cleanup SDL resources.
     */
    ~Renderer();

    /**
     * @brief Clear the framebuffer with a color.
     */
    void clear(const Color& color = Color::Black());

    /**
     * @brief Render a mesh to the framebuffer.
     * @param mesh The mesh to render (in world coordinates).
     * @param line_color Color for rendering (currently unused in filled rendering).
     */
    void renderMesh(const Mesh& mesh, const Color& line_color = Color::White());

    /**
     * @brief Render a shape at a specific position and orientation.
     * @param mesh The shape's mesh in local coordinates.
     * @param transform The transform specifying position and orientation.
     * @param color Color for rendering the shape.
     */
    void renderShape(const Mesh& mesh, const Transform& transform, const Color& color = Color::White());

    /**
     * @brief Render a rigid body.
     * @param body The rigid body to render.
     */
    void renderBody(const RigidBody& body);

    /**
     * @brief Draw coordinate axes at origin.
     * @param length Length of each axis.
     */
    void renderWorldFrame(real_t length = 1.0);

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