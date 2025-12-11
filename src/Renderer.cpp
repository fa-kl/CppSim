
#include "Renderer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace sim
{

#pragma region FrameBuffer

FrameBuffer::FrameBuffer(int_t width, int_t height)
    : m_width(width), m_height(height), m_pixel_colors(static_cast<uint_t>(width * height))
{
}

bool FrameBuffer::isInFrame(const PixelPosition& px) const
{
  return 0 <= px.x && px.x < m_width && 0 <= px.y && px.y < m_height;
}

Color& FrameBuffer::operator[](PixelPosition px)
{
  if (isInFrame(px)) {
    return m_pixel_colors[static_cast<uint_t>(px.y * m_width + px.x)];
  }
  std::string pixel = "(" + std::to_string(px.x) + ", " + std::to_string(px.y) + ")";
  std::string window_size = "(" + std::to_string(m_width) + ", " + std::to_string(m_height) + ")";
  std::string msg = "Pixel " + pixel + "out of range for window of size " + window_size + ".";
  throw std::out_of_range(msg);
}

const Color& FrameBuffer::operator[](PixelPosition px) const
{
  if (isInFrame(px)) {
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

#pragma endregion

#pragma region MeshTransform

Mesh operator*(const Transform<3>& transform, const Mesh& mesh)
{
  Mesh result = mesh;
  for (Triangle& triangle : result) {
    for (Vertex& vert : triangle) {
      vert.point = transform * vert.point;
    }
  }
  return result;
}

Mesh cullMesh(const Mesh& input_mesh, const Camera& camera)
{
  Mesh result;
  result.reserve(input_mesh.size());
  for (const auto& triangle : input_mesh) {
    bool any_in_view = false;
    for (const auto& vertex : triangle) {
      if (camera.isInView(vertex.point)) {
        any_in_view = true;
        break;
      }
    }
    if (any_in_view) {
      result.push_back(triangle);
    }
  }
  return result;
}

#pragma endregion

#pragma region Renderer

Renderer::Renderer(int_t width, int_t height, const std::string& title)
    : m_framebuffer(width, height), m_width(width), m_height(height), m_should_close(false)
{
  if (SDL_Init(SDL_INIT_VIDEO) < 0) {
    throw std::runtime_error(std::string("SDL could not initialize! SDL_Error: ") + SDL_GetError());
  }
  m_window = SDL_CreateWindow(title.c_str(),
                              SDL_WINDOWPOS_CENTERED,
                              SDL_WINDOWPOS_CENTERED,
                              static_cast<int>(width),
                              static_cast<int>(height),
                              SDL_WINDOW_SHOWN);
  if (!m_window) {
    SDL_Quit();
    throw std::runtime_error(std::string("Window could not be created! SDL_Error: ") + SDL_GetError());
  }
  m_sdl_renderer = SDL_CreateRenderer(m_window, -1, SDL_RENDERER_ACCELERATED);
  if (!m_sdl_renderer) {
    SDL_DestroyWindow(m_window);
    SDL_Quit();
    throw std::runtime_error(std::string("Renderer could not be created! SDL_Error: ") + SDL_GetError());
  }
  m_texture = SDL_CreateTexture(m_sdl_renderer,
                                SDL_PIXELFORMAT_RGBA32,
                                SDL_TEXTUREACCESS_STREAMING,
                                static_cast<int>(width),
                                static_cast<int>(height));
  if (!m_texture) {
    SDL_DestroyRenderer(m_sdl_renderer);
    SDL_DestroyWindow(m_window);
    SDL_Quit();
    throw std::runtime_error(std::string("Texture could not be created! SDL_Error: ") + SDL_GetError());
  }
}

Renderer::~Renderer()
{
  if (m_texture)
    SDL_DestroyTexture(m_texture);
  if (m_sdl_renderer)
    SDL_DestroyRenderer(m_sdl_renderer);
  if (m_window)
    SDL_DestroyWindow(m_window);
  SDL_Quit();
}

void Renderer::clear(const Color& color)
{
  m_framebuffer.clear(color);
}

void Renderer::drawLine(int x0, int y0, int x1, int y1, const Color& color)
{
  int dx = std::abs(x1 - x0), dy = std::abs(y1 - y0);
  int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1, err = dx - dy;
  while (true) {
    if (x0 >= 0 && x0 < m_width && y0 >= 0 && y0 < m_height) {
      m_framebuffer[{x0, y0}] = color;
    }
    if (x0 == x1 && y0 == y1) {
      break;
    }
    int e2 = 2 * err;
    if (e2 > -dy) {
      err -= dy;
      x0 += sx;
    }
    if (e2 < dx) {
      err += dx;
      y0 += sy;
    }
  }
}

void Renderer::fillTriangle(int x0, int y0, int x1, int y1, int x2, int y2, const Color& color)
{
  if (y0 > y1) {
    std::swap(y0, y1);
    std::swap(x0, x1);
  }
  if (y0 > y2) {
    std::swap(y0, y2);
    std::swap(x0, x2);
  }
  if (y1 > y2) {
    std::swap(y1, y2);
    std::swap(x1, x2);
  }
  auto interp = [](int y, int y0, int x0, int y1, int x1) {
    return y1 == y0 ? x0 : x0 + (x1 - x0) * (y - y0) / (y1 - y0);
  };
  for (int y = y0; y <= y2; ++y) {
    int xa = interp(y, y0, x0, y2, x2), xb = y <= y1 ? interp(y, y0, x0, y1, x1) : interp(y, y1, x1, y2, x2);
    if (xa > xb) {
      std::swap(xa, xb);
    }
    for (int x = xa; x <= xb; ++x) {
      if (x >= 0 && x < m_width && y >= 0 && y < m_height) {
        m_framebuffer[{x, y}] = color;
      }
    }
  }
}

std::pair<int, int> Renderer::ndcToScreen(const Vector& ndc) const
{
  return {static_cast<int>((ndc[0] + 1.0) * 0.5 * m_width), static_cast<int>((1.0 - ndc[1]) * 0.5 * m_height)};
}

void Renderer::renderMesh(const Mesh& mesh, const Camera& camera, const Color& line_color)
{
  // First pass: fill triangles with vertex colors
  for (const auto& tri : mesh) {
    try {
      auto [x0, y0] = ndcToScreen(camera.projectWorldToNDC(tri[0].point));
      auto [x1, y1] = ndcToScreen(camera.projectWorldToNDC(tri[1].point));
      auto [x2, y2] = ndcToScreen(camera.projectWorldToNDC(tri[2].point));
      // Use the color from the first vertex (or could average all three)
      fillTriangle(x0, y0, x1, y1, x2, y2, tri[0].color);
    } catch (...) {
    }
  }
  // Second pass: draw wireframe lines
  for (const auto& tri : mesh) {
    try {
      auto [x0, y0] = ndcToScreen(camera.projectWorldToNDC(tri[0].point));
      auto [x1, y1] = ndcToScreen(camera.projectWorldToNDC(tri[1].point));
      auto [x2, y2] = ndcToScreen(camera.projectWorldToNDC(tri[2].point));
      drawLine(x0, y0, x1, y1, line_color);
      drawLine(x1, y1, x2, y2, line_color);
      drawLine(x2, y2, x0, y0, line_color);
    } catch (...) {
    }
  }
}

void Renderer::drawAxes(const Camera& camera, real_t length)
{
  Vector o = {0, 0, 0}, x = {length, 0, 0}, y = {0, length, 0}, z = {0, 0, length};
  try {
    if (camera.isInView(o) || camera.isInView(x)) {
      auto [ox, oy] = ndcToScreen(camera.projectWorldToNDC(o));
      auto [xx, xy] = ndcToScreen(camera.projectWorldToNDC(x));
      drawLine(ox, oy, xx, xy, RED);
    }
    if (camera.isInView(o) || camera.isInView(y)) {
      auto [ox, oy] = ndcToScreen(camera.projectWorldToNDC(o));
      auto [yx, yy] = ndcToScreen(camera.projectWorldToNDC(y));
      drawLine(ox, oy, yx, yy, GREEN);
    }
    if (camera.isInView(o) || camera.isInView(z)) {
      auto [ox, oy] = ndcToScreen(camera.projectWorldToNDC(o));
      auto [zx, zy] = ndcToScreen(camera.projectWorldToNDC(z));
      drawLine(ox, oy, zx, zy, BLUE);
    }
  } catch (...) {
  }
}

void Renderer::present()
{
  SDL_UpdateTexture(m_texture, nullptr, m_framebuffer.getRawData(), static_cast<int>(m_width) * 4);
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
  return m_width;
}
int_t Renderer::getHeight() const
{
  return m_height;
}
bool Renderer::shouldClose()
{
  return m_should_close;
}

void Renderer::pollEvents()
{
  SDL_Event e;
  while (SDL_PollEvent(&e) != 0) {
    if (e.type == SDL_QUIT || (e.type == SDL_KEYDOWN && e.key.keysym.sym == SDLK_ESCAPE))
      m_should_close = true;
  }
}

#pragma endregion

}  // namespace sim
