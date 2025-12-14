
#include "Renderer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace sim
{

Renderer::Renderer(const Camera& camera, const std::string& title)
    : m_framebuffer(camera.getWidth(), camera.getHeight()), m_camera(camera), m_should_close(false)
{
    if (SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        throw std::runtime_error(std::string("SDL could not initialize! SDL_Error: ") + SDL_GetError());
    }
    m_window = SDL_CreateWindow(title.c_str(),
                                SDL_WINDOWPOS_CENTERED,
                                SDL_WINDOWPOS_CENTERED,
                                static_cast<int>(camera.getWidth()),
                                static_cast<int>(camera.getHeight()),
                                SDL_WINDOW_SHOWN);
    if (!m_window)
    {
        SDL_Quit();
        throw std::runtime_error(std::string("Window could not be created! SDL_Error: ") + SDL_GetError());
    }
    m_sdl_renderer = SDL_CreateRenderer(m_window, -1, SDL_RENDERER_ACCELERATED);
    if (!m_sdl_renderer)
    {
        SDL_DestroyWindow(m_window);
        SDL_Quit();
        throw std::runtime_error(std::string("Renderer could not be created! SDL_Error: ") + SDL_GetError());
    }
    m_texture = SDL_CreateTexture(m_sdl_renderer,
                                  SDL_PIXELFORMAT_RGBA32,
                                  SDL_TEXTUREACCESS_STREAMING,
                                  static_cast<int>(camera.getWidth()),
                                  static_cast<int>(camera.getHeight()));
    if (!m_texture)
    {
        SDL_DestroyRenderer(m_sdl_renderer);
        SDL_DestroyWindow(m_window);
        SDL_Quit();
        throw std::runtime_error(std::string("Texture could not be created! SDL_Error: ") + SDL_GetError());
    }
}

Renderer::~Renderer()
{
    if (m_texture)
    {
        SDL_DestroyTexture(m_texture);
    }
    if (m_sdl_renderer)
    {
        SDL_DestroyRenderer(m_sdl_renderer);
    }
    if (m_window)
    {
        SDL_DestroyWindow(m_window);
    }
    SDL_Quit();
}

void Renderer::drawLine(const Vector& v1, const Vector& v2, const Color& color)
{
    // Project world points to NDC
    Vector ndc1 = m_camera.projectWorldToNDC(v1);
    Vector ndc2 = m_camera.projectWorldToNDC(v2);

    // Convert NDC to screen coordinates
    // NDC: x,y in [-1, 1] -> Screen: x in [0, width], y in [0, height]
    // Note: y-axis is flipped (NDC +y is up, screen +y is down)
    Vector screen1 = {(ndc1[0] + 1.0) * 0.5 * static_cast<real_t>(m_camera.getWidth()),
                      (1.0 - ndc1[1]) * 0.5 * static_cast<real_t>(m_camera.getHeight())};
    Vector screen2 = {(ndc2[0] + 1.0) * 0.5 * static_cast<real_t>(m_camera.getWidth()),
                      (1.0 - ndc2[1]) * 0.5 * static_cast<real_t>(m_camera.getHeight())};

    // Bresenham's line algorithm using real-valued vectors for accuracy
    Vector delta = screen2 - screen1;
    Vector abs_delta = abs(delta);

    // Determine the number of steps
    int_t steps = static_cast<int_t>(max(abs_delta)) + 1;

    if (steps <= 1)
    {
        // Single pixel
        PixelCoordinate px(screen1);
        // Check bounds manually
        if (px.x >= 0 && px.x < m_camera.getWidth() && px.y >= 0 && px.y < m_camera.getHeight())
        {
            m_framebuffer[px] = color;
        }
        return;
    }

    // Interpolate along the line
    for (int_t i = 0; i <= steps; ++i)
    {
        real_t t = static_cast<real_t>(i) / static_cast<real_t>(steps);
        Vector current = screen1 + t * delta;
        PixelCoordinate px(current);

        // Check bounds manually
        if (px.x >= 0 && px.x < m_camera.getWidth() && px.y >= 0 && px.y < m_camera.getHeight())
        {
            m_framebuffer[px] = color;
        }
    }
}

void Renderer::fillTriangle(const Vertex& v1, const Vertex& v2, const Vertex& v3)
{
    // Project vertices to NDC
    Vector ndc1 = m_camera.projectWorldToNDC(v1.position);
    Vector ndc2 = m_camera.projectWorldToNDC(v2.position);
    Vector ndc3 = m_camera.projectWorldToNDC(v3.position);

    // Extract depth values from NDC z-coordinates (in range [-1, 1])
    real_t depth1 = ndc1[2];
    real_t depth2 = ndc2[2];
    real_t depth3 = ndc3[2];

    // Convert NDC to screen coordinates (keeping as Vectors for precision)
    Vector screen1 = {(ndc1[0] + 1.0) * 0.5 * static_cast<real_t>(m_camera.getWidth()),
                      (1.0 - ndc1[1]) * 0.5 * static_cast<real_t>(m_camera.getHeight())};
    Vector screen2 = {(ndc2[0] + 1.0) * 0.5 * static_cast<real_t>(m_camera.getWidth()),
                      (1.0 - ndc2[1]) * 0.5 * static_cast<real_t>(m_camera.getHeight())};
    Vector screen3 = {(ndc3[0] + 1.0) * 0.5 * static_cast<real_t>(m_camera.getWidth()),
                      (1.0 - ndc3[1]) * 0.5 * static_cast<real_t>(m_camera.getHeight())};

    // Find bounding box
    real_t min_x = std::min({screen1[0], screen2[0], screen3[0]});
    real_t max_x = std::max({screen1[0], screen2[0], screen3[0]});
    real_t min_y = std::min({screen1[1], screen2[1], screen3[1]});
    real_t max_y = std::max({screen1[1], screen2[1], screen3[1]});

    // Clamp to screen bounds
    int_t start_x = std::max(static_cast<int_t>(0), static_cast<int_t>(std::floor(min_x)));
    int_t end_x = std::min(static_cast<int_t>(m_camera.getWidth() - 1), static_cast<int_t>(std::ceil(max_x)));
    int_t start_y = std::max(static_cast<int_t>(0), static_cast<int_t>(std::floor(min_y)));
    int_t end_y = std::min(static_cast<int_t>(m_camera.getHeight() - 1), static_cast<int_t>(std::ceil(max_y)));

    // Precompute edge vectors for barycentric coordinates
    Vector edge0 = screen2 - screen1;
    Vector edge1 = screen3 - screen1;

    real_t dot00 = dot(edge0, edge0);
    real_t dot01 = dot(edge0, edge1);
    real_t dot11 = dot(edge1, edge1);
    real_t inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01);

    // Scan through bounding box
    for (int_t py = start_y; py <= end_y; ++py)
    {
        for (int_t px = start_x; px <= end_x; ++px)
        {
            // Use real coordinates for the pixel center
            Vector p = {static_cast<real_t>(px) + 0.5, static_cast<real_t>(py) + 0.5};
            Vector edge2 = p - screen1;

            // Compute barycentric coordinates
            real_t dot02 = dot(edge0, edge2);
            real_t dot12 = dot(edge1, edge2);
            real_t bary_u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
            real_t bary_v = (dot00 * dot12 - dot01 * dot02) * inv_denom;
            real_t bary_w = 1.0 - bary_u - bary_v;

            // Check if point is inside triangle
            if (bary_u >= 0.0 && bary_v >= 0.0 && bary_w >= 0.0)
            {
                // Interpolate depth using barycentric coordinates
                real_t pixel_depth = bary_w * depth1 + bary_u * depth2 + bary_v * depth3;

                // Depth test: only draw if this pixel is closer than what's currently in the buffer
                if (m_framebuffer.testAndSetDepth({px, py}, pixel_depth))
                {
                    m_framebuffer[{px, py}] = interpolate(v1.color, v2.color, v3.color, {bary_w, bary_u, bary_v});
                }
            }
        }
    }
}

void Renderer::clear(const Color& color)
{
    m_framebuffer.clear(color);
}

void Renderer::renderMesh(const Mesh& mesh, const Color& line_color)
{
    const auto& vertices = mesh.getVertices();
    const auto& faces = mesh.getFaceIndices();

    for (const auto& face : faces)
    {
        try
        {
            // Fill the triangle
            fillTriangle(vertices[face[0]], vertices[face[1]], vertices[face[2]]);

            // Draw the edges
            drawLine(vertices[face[0]].position, vertices[face[1]].position, line_color);
            drawLine(vertices[face[1]].position, vertices[face[2]].position, line_color);
            drawLine(vertices[face[2]].position, vertices[face[0]].position, line_color);
        }
        catch (...)
        {
            // Skip triangles that can't be rendered
        }
    }
}

void Renderer::renderShape(const Mesh& mesh, const Transform& transform, const Color& color)
{
    // Transform mesh to world coordinates
    Mesh world_mesh = transform * mesh;

    // Set all vertex colors to the specified color
    world_mesh.setUniformColor(color);

    // Cull mesh based on camera frustum
    Mesh culled_mesh = cullMesh(world_mesh, m_camera);

    // Render the culled mesh
    renderMesh(culled_mesh, color);
}

void Renderer::renderWorldFrame(real_t length)
{
    Vector o = {0, 0, 0}, x = {length, 0, 0}, y = {0, length, 0}, z = {0, 0, length};
    try
    {
        const Vector origin = m_camera.projectWorldToCamera({0, 0, 0});
        const Vector x_axis = m_camera.projectWorldToCamera({1, 0, 0});
        const Vector y_axis = m_camera.projectWorldToCamera({0, 1, 0});
        const Vector z_axis = m_camera.projectWorldToCamera({0, 0, 1});
        if (m_camera.isInView(o) || m_camera.isInView(x))
        {
            drawLine(origin, x_axis, Color::Red());
        }
        if (m_camera.isInView(o) || m_camera.isInView(y))
        {
            drawLine(origin, y_axis, Color::Green());
        }
        if (m_camera.isInView(o) || m_camera.isInView(z))
        {
            drawLine(origin, z_axis, Color::Blue());
        }
    }
    catch (...)
    {
    }
}

void Renderer::present()
{
    SDL_UpdateTexture(m_texture, nullptr, m_framebuffer.getRawData(), static_cast<int>(m_camera.getWidth()) * 4);
    SDL_RenderClear(m_sdl_renderer);
    SDL_RenderCopy(m_sdl_renderer, m_texture, nullptr, nullptr);
    SDL_RenderPresent(m_sdl_renderer);
}

FrameBuffer& Renderer::getFrameBuffer()
{
    return m_framebuffer;
}

int_t Renderer::getWidth() const
{
    return m_camera.getWidth();
}

int_t Renderer::getHeight() const
{
    return m_camera.getHeight();
}

bool Renderer::shouldClose()
{
    return m_should_close;
}

void Renderer::pollEvents()
{
    SDL_Event e;
    while (SDL_PollEvent(&e) != 0)
    {
        if (e.type == SDL_QUIT || (e.type == SDL_KEYDOWN && e.key.keysym.sym == SDLK_ESCAPE))
            m_should_close = true;
    }
}

void Renderer::renderBody(const RigidBody& body)
{
    renderShape(body.getMesh(), body.getTransform(), body.getMaterial().color);
}

}  // namespace sim
