
#include "Renderer.hpp"

#include <stdexcept>
#include <string>

namespace sim
{

#pragma region PixelPosition

PixelPosition operator+(PixelPosition px1, const PixelPosition& px2)
{
  px1.x += px2.x;
  px1.y += px2.y;
  return px1;
}

PixelPosition operator-(PixelPosition px1, const PixelPosition& px2)
{
  px1.x -= px2.x;
  px1.y -= px2.y;
  return px1;
}

PixelPosition operator*(PixelPosition px, const int_t& k)
{
  px.x *= static_cast<uint_t>(k);
  px.y *= static_cast<uint_t>(k);
  return px;
}

PixelPosition operator*(const int_t& k, PixelPosition px)
{
  px.x *= static_cast<uint_t>(k);
  px.y *= static_cast<uint_t>(k);
  return px;
}

#pragma endregion

#pragma region Color

Color color(real_t r, real_t g, real_t b, real_t opacity)
{
  if (0 < r || r < 1) {
    throw std::invalid_argument("All values must be within [0; 1].");
  }
  if (0 < g || g < 1) {
    throw std::invalid_argument("All values must be within [0; 1].");
  }
  if (0 < b || b < 1) {
    throw std::invalid_argument("All values must be within [0; 1].");
  }
  if (0 < opacity || opacity < 1) {
    throw std::invalid_argument("All values must be within [0; 1].");
  }
  uint8_t r8 = static_cast<uint8_t>(round(255.0 * r));
  uint8_t g8 = static_cast<uint8_t>(round(255.0 * g));
  uint8_t b8 = static_cast<uint8_t>(round(255.0 * b));
  uint8_t a8 = static_cast<uint8_t>(round(255.0 * opacity));
  return {r8, g8, b8, a8};
}

const Color white = {255, 255, 255, 255};
const Color black = {0, 0, 0, 255};
const Color transparent = {0, 0, 0, 0};
const Color red = {255, 0, 0, 255};
const Color green = {0, 255, 0, 255};
const Color blue = {0, 0, 255, 255};
const Color yellow = {255, 255, 0, 255};
const Color cyan = {0, 255, 255, 255};
const Color magenta = {255, 0, 255, 255};
const Color orange = {255, 165, 0, 255};
const Color purple = {128, 0, 128, 255};
const Color pink = {255, 192, 203, 255};
const Color lightRed = {255, 128, 128, 255};
const Color lightGreen = {144, 238, 144, 255};
const Color lightBlue = {173, 216, 230, 255};
const Color lightYellow = {255, 255, 224, 255};
const Color lightCyan = {224, 255, 255, 255};
const Color lightMagenta = {255, 128, 255, 255};
const Color lightOrange = {255, 200, 124, 255};
const Color lightPurple = {200, 162, 200, 255};
const Color lightPink = {255, 220, 225, 255};
const Color lightGray = {211, 211, 211, 255};
const Color darkRed = {139, 0, 0, 255};
const Color darkGreen = {0, 100, 0, 255};
const Color darkBlue = {0, 0, 139, 255};
const Color darkYellow = {204, 204, 0, 255};
const Color darkCyan = {0, 139, 139, 255};
const Color darkMagenta = {139, 0, 139, 255};
const Color darkOrange = {255, 140, 0, 255};
const Color darkPurple = {75, 0, 130, 255};
const Color darkPink = {231, 84, 128, 255};
const Color darkGray = {169, 169, 169, 255};
const Color gray = {128, 128, 128, 255};
const Color silver = {192, 192, 192, 255};
const Color dimGray = {105, 105, 105, 255};
const Color redAlpha = {255, 0, 0, 127};
const Color greenAlpha = {0, 255, 0, 127};
const Color blueAlpha = {0, 0, 255, 127};
const Color yellowAlpha = {255, 255, 0, 127};
const Color cyanAlpha = {0, 255, 255, 127};
const Color magentaAlpha = {255, 0, 255, 127};
const Color whiteAlpha = {255, 255, 255, 127};
const Color blackAlpha = {0, 0, 0, 127};

#pragma endregion

#pragma region FrameBuffer

FrameBuffer::FrameBuffer(uint_t width, uint_t height) : m_width(width), m_height(height), m_pixel_colors(width * height)
{
}

Color& FrameBuffer::operator[](PixelPosition px)
{
  if (px.x < m_width && px.y < m_height) {
    return m_pixel_colors[px.y * m_width + px.x];
  }
  std::string pixel = "(" + std::to_string(px.x) + ", " + std::to_string(px.y) + ")";
  std::string window_size = "(" + std::to_string(m_width) + ", " + std::to_string(m_height) + ")";
  std::string msg = "Pixel " + pixel + "out of range for window of size " + window_size + ".";
  throw std::out_of_range(msg);
}

const Color& FrameBuffer::operator[](PixelPosition px) const
{
  if (px.x < m_width && px.y < m_height) {
    return m_pixel_colors[px.y * m_width + px.x];
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

uint_t FrameBuffer::getWidth() const
{
  return m_width;
}
uint_t FrameBuffer::getHeight() const
{
  return m_height;
}

const void* FrameBuffer::getRawData() const
{
  return m_pixel_colors.data();
}

#pragma endregion

#pragma region Renderer

Renderer::Renderer(FrameBuffer& fb) : m_buffer(fb) {}

void Renderer::drawLine(const PixelPosition& start, const PixelPosition& end, uint_t width, const Color& color)
{
  int_t x0 = static_cast<int_t>(start.x);
  int_t y0 = static_cast<int_t>(start.y);
  int_t x1 = static_cast<int_t>(end.x);
  int_t y1 = static_cast<int_t>(end.y);
  int_t dx = std::abs(x1 - x0);
  int_t dy = std::abs(y1 - y0);
  int_t sx = (x0 < x1) ? 1 : -1;
  int_t sy = (y0 < y1) ? 1 : -1;
  int_t err = dx - dy;

  width = 2 * (width - 1) + 1;
  int_t halfWidth = static_cast<int_t>(width) / 2;

  while (true) {
    for (int_t w = -halfWidth; w <= halfWidth; w++) {
      for (int_t h = -halfWidth; h <= halfWidth; h++) {
        int_t px = x0 + w;
        int_t py = y0 + h;
        if (px >= 0 && py >= 0 && static_cast<uint_t>(px) < m_buffer.getWidth() &&
            static_cast<uint_t>(py) < m_buffer.getHeight()) {
          PixelPosition pixelpos = {static_cast<uint_t>(px), static_cast<uint_t>(py)};
          m_buffer[pixelpos] = color;
        }
      }
    }
    if (x0 == x1 && y0 == y1) {
      break;
    }
    int_t e2 = 2 * err;
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

void Renderer::drawCircle(const PixelPosition& center, const real_t& radius, const Color& color)
{
  int_t cx = static_cast<int_t>(center.x);
  int_t cy = static_cast<int_t>(center.y);
  int_t r = static_cast<int_t>(radius);
  for (int_t y = -r; y <= r; y++) {
    for (int_t x = -r; x <= r; x++) {
      if (x * x + y * y <= r * r) {
        int_t px = cx + x;
        int_t py = cy + y;
        if (px >= 0 && py >= 0 && static_cast<uint_t>(px) < m_buffer.getWidth() &&
            static_cast<uint_t>(py) < m_buffer.getHeight()) {
          PixelPosition pixelpos = {static_cast<uint_t>(px), static_cast<uint_t>(py)};
          m_buffer[pixelpos] = color;
        }
      }
    }
  }
}

void Renderer::drawPolygon(const std::vector<PixelPosition>& vertices, const uint_t& lineWidth, const Color& color)
{
  if (vertices.size() < 3) {
    return;
  }
  for (size_t i = 0; i < vertices.size(); i++) {
    size_t next = (i + 1) % vertices.size();
    drawLine(vertices[i], vertices[next], lineWidth, color);
  }
}

#pragma endregion

#pragma region Window

Window::Window(uint_t w, uint_t h) : width(w), height(h)
{
  if (SDL_Init(SDL_INIT_VIDEO) < 0) {
    throw std::runtime_error("SDL init failed: " + std::string(SDL_GetError()));
  }
  window = SDL_CreateWindow("Custom 2D Renderer",
                            SDL_WINDOWPOS_CENTERED,
                            SDL_WINDOWPOS_CENTERED,
                            static_cast<int>(width),
                            static_cast<int>(height),
                            SDL_WINDOW_SHOWN);

  if (!window) {
    throw std::runtime_error("Window creation failed: " + std::string(SDL_GetError()));
  }
  sdlRenderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
  if (!sdlRenderer) {
    throw std::runtime_error("Renderer creation failed: " + std::string(SDL_GetError()));
  }
  texture = SDL_CreateTexture(sdlRenderer,
                              SDL_PIXELFORMAT_ABGR8888,
                              SDL_TEXTUREACCESS_STREAMING,
                              static_cast<int>(width),
                              static_cast<int>(height));
  if (!texture) {
    throw std::runtime_error("Texture creation failed: " + std::string(SDL_GetError()));
  }
}

Window::~Window()
{
  if (texture)
    SDL_DestroyTexture(texture);
  if (sdlRenderer)
    SDL_DestroyRenderer(sdlRenderer);
  if (window)
    SDL_DestroyWindow(window);
  SDL_Quit();
}

void Window::display(const FrameBuffer& frameBuffer)
{
  SDL_UpdateTexture(texture, nullptr, frameBuffer.getRawData(), static_cast<int>(width * 4));
  SDL_SetRenderDrawColor(sdlRenderer, 0, 0, 0, 255);
  SDL_RenderClear(sdlRenderer);
  SDL_RenderCopy(sdlRenderer, texture, nullptr, nullptr);
  SDL_RenderPresent(sdlRenderer);
}

bool Window::shouldClose()
{
  SDL_Event event;
  while (SDL_PollEvent(&event)) {
    if (event.type == SDL_QUIT || (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_ESCAPE)) {
      return true;
    }
  }
  return false;
}

#pragma endregion

}  // namespace sim