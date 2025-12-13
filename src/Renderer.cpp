
#include "Renderer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace sim
{

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
  if (m_texture) {
    SDL_DestroyTexture(m_texture);
  }
  if (m_sdl_renderer) {
    SDL_DestroyRenderer(m_sdl_renderer);
  }
  if (m_window) {
    SDL_DestroyWindow(m_window);
  }
  SDL_Quit();
}

void Renderer::clear(const Color& color)
{
  m_framebuffer.clear(color);
}

void Renderer::drawLine(PixelCoordinate px1, PixelCoordinate px2, const Color& color)
{
  int_t dx = std::abs(px2.x - px1.x), dy = std::abs(px2.y - px1.y);
  int_t sx = px1.x < px2.x ? 1 : -1, sy = px1.y < px2.y ? 1 : -1, err = dx - dy;
  while (true) {
    if (px1.x >= 0 && px1.x < m_width && px1.y >= 0 && px1.y < m_height) {
      m_framebuffer[{px1.x, px1.y}] = color;
    }
    if (px1.x == px2.x && px1.y == px2.y) {
      break;
    }
    int_t e2 = 2 * err;
    if (e2 > -dy) {
      err -= dy;
      px1.x += sx;
    }
    if (e2 < dx) {
      err += dx;
      px1.y += sy;
    }
  }
}

void Renderer::fillTriangle(PixelCoordinate px1, PixelCoordinate px2, PixelCoordinate px3, const Color& color)
{
  if (px1.y > px2.y) {
    std::swap(px1.y, px2.y);
    std::swap(px1.x, px2.x);
  }
  if (px1.y > px3.y) {
    std::swap(px1.y, px3.y);
    std::swap(px1.x, px3.x);
  }
  if (px2.y > px3.y) {
    std::swap(px2.y, px3.y);
    std::swap(px2.x, px3.x);
  }
  auto interp = [](int_t y, const PixelCoordinate& p1, const PixelCoordinate& p2) {
    return p2.y == p1.y ? p1.x : p1.x + (p2.x - p1.x) * (y - p1.y) / (p2.y - p1.y);
  };
  for (int_t y = px1.y; y <= px3.y; ++y) {
    int_t xa = interp(y, px1, px3), xb = y <= px2.y ? interp(y, px1, px2) : interp(y, px2, px3);
    if (xa > xb) {
      std::swap(xa, xb);
    }
    for (int_t x = xa; x <= xb; ++x) {
      if (x >= 0 && x < m_width && y >= 0 && y < m_height) {
        m_framebuffer[{x, y}] = color;
      }
    }
  }
}

PixelCoordinate Renderer::ndcToScreen(const Vector& ndc) const
{
  return {static_cast<int_t>((ndc[0] + 1.0) * 0.5 * static_cast<real_t>(m_width)),
          static_cast<int_t>((1.0 - ndc[1]) * 0.5 * static_cast<real_t>(m_height))};
}

void Renderer::renderMesh(const Mesh& mesh, const Camera& camera, const Color& line_color)
{
  for (const Triangle& tri : mesh) {
    try {
      const PixelCoordinate px1 = ndcToScreen(camera.projectWorldToNDC(tri[0].point));
      const PixelCoordinate px2 = ndcToScreen(camera.projectWorldToNDC(tri[1].point));
      const PixelCoordinate px3 = ndcToScreen(camera.projectWorldToNDC(tri[2].point));
      fillTriangle(px1, px2, px3, mean({tri[0].color, tri[1].color, tri[2].color}));
    } catch (...) {
    }
  }
  return;
  for (const Triangle& tri : mesh) {
    try {
      const PixelCoordinate px1 = ndcToScreen(camera.projectWorldToNDC(tri[0].point));
      const PixelCoordinate px2 = ndcToScreen(camera.projectWorldToNDC(tri[1].point));
      const PixelCoordinate px3 = ndcToScreen(camera.projectWorldToNDC(tri[2].point));
      drawLine(px1, px2, line_color);
      drawLine(px2, px3, line_color);
      drawLine(px3, px1, line_color);
    } catch (...) {
    }
  }
}

void Renderer::drawAxes(const Camera& camera, real_t length)
{
  Vector o = {0, 0, 0}, x = {length, 0, 0}, y = {0, length, 0}, z = {0, 0, length};
  try {
    const PixelCoordinate origin = ndcToScreen(camera.projectWorldToNDC(o));
    if (camera.isInView(o) || camera.isInView(x)) {
      const PixelCoordinate x_axis = ndcToScreen(camera.projectWorldToNDC(x));
      drawLine(origin, x_axis, RED);
    }
    if (camera.isInView(o) || camera.isInView(y)) {
      const PixelCoordinate y_axis = ndcToScreen(camera.projectWorldToNDC(y));
      drawLine(origin, y_axis, GREEN);
    }
    if (camera.isInView(o) || camera.isInView(z)) {
      const PixelCoordinate z_axis = ndcToScreen(camera.projectWorldToNDC(z));
      drawLine(origin, z_axis, BLUE);
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

}  // namespace sim
