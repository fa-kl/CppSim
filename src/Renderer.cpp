
#include "Renderer.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace sim
{

#pragma region Color

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

FrameBuffer::FrameBuffer(int_t width, int_t height)
    : m_width(width), m_height(height), m_pixel_colors(static_cast<uint_t>(width * height))
{
}

Color& FrameBuffer::operator[](PixelPosition px)
{
  if (0 <= px.x && px.x < m_width && 0 <= px.y && px.y < m_height) {
    return m_pixel_colors[static_cast<uint_t>(px.y * m_width + px.x)];
  }
  std::string pixel = "(" + std::to_string(px.x) + ", " + std::to_string(px.y) + ")";
  std::string window_size = "(" + std::to_string(m_width) + ", " + std::to_string(m_height) + ")";
  std::string msg = "Pixel " + pixel + "out of range for window of size " + window_size + ".";
  throw std::out_of_range(msg);
}

const Color& FrameBuffer::operator[](PixelPosition px) const
{
  if (0 <= px.x && px.x < m_width && 0 <= px.y && px.y < m_height) {
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

#pragma region Renderer

PixelPosition Renderer::worldToPixel(const Vector& worldPos) const
{
  if (worldPos.length() != 2) {
    throw std::invalid_argument("Vector must be 2D for rendering");
  }
  int_t px = m_center.x + static_cast<int_t>(std::round(worldPos[0] * m_scale));
  int_t py = m_center.y - static_cast<int_t>(std::round(worldPos[1] * m_scale));
  return {px, py};
}

Renderer::Renderer(int_t width, int_t height, real_t scale)
    : m_buffer(width, height),
      m_scale(scale),
      m_center({width / 2, height / 2}),
      m_grid_spacing(static_cast<int_t>(scale)),
      m_width(width),
      m_height(height)
{
  if (SDL_Init(SDL_INIT_VIDEO) < 0) {
    throw std::runtime_error("SDL init failed: " + std::string(SDL_GetError()));
  }
  m_window = SDL_CreateWindow("CppSim Renderer",
                              SDL_WINDOWPOS_CENTERED,
                              SDL_WINDOWPOS_CENTERED,
                              static_cast<int>(width),
                              static_cast<int>(height),
                              SDL_WINDOW_SHOWN);

  if (!m_window) {
    throw std::runtime_error("Window creation failed: " + std::string(SDL_GetError()));
  }
  m_SDLRenderer = SDL_CreateRenderer(m_window, -1, SDL_RENDERER_ACCELERATED);
  if (!m_SDLRenderer) {
    throw std::runtime_error("Renderer creation failed: " + std::string(SDL_GetError()));
  }
  m_texture = SDL_CreateTexture(m_SDLRenderer,
                                SDL_PIXELFORMAT_ABGR8888,
                                SDL_TEXTUREACCESS_STREAMING,
                                static_cast<int>(width),
                                static_cast<int>(height));
  if (!m_texture) {
    throw std::runtime_error("Texture creation failed: " + std::string(SDL_GetError()));
  }
}

Renderer::~Renderer()
{
  if (m_texture)
    SDL_DestroyTexture(m_texture);
  if (m_SDLRenderer)
    SDL_DestroyRenderer(m_SDLRenderer);
  if (m_window)
    SDL_DestroyWindow(m_window);
  SDL_Quit();
}

void Renderer::drawLine(const Vector& start, const Vector& end, uint_t widthPixels, const Color& color)
{
  PixelPosition startPx = worldToPixel(start);
  PixelPosition endPx = worldToPixel(end);

  int_t x0 = startPx.x;
  int_t y0 = startPx.y;
  int_t x1 = endPx.x;
  int_t y1 = endPx.y;
  int_t dx = std::abs(x1 - x0);
  int_t dy = std::abs(y1 - y0);
  int_t sx = (x0 < x1) ? 1 : -1;
  int_t sy = (y0 < y1) ? 1 : -1;
  int_t err = dx - dy;

  widthPixels = 2 * (widthPixels - 1) + 1;
  int_t halfWidth = static_cast<int_t>(widthPixels) / 2;

  while (true) {
    for (int_t w = -halfWidth; w <= halfWidth; w++) {
      for (int_t h = -halfWidth; h <= halfWidth; h++) {
        int_t px = x0 + w;
        int_t py = y0 + h;
        if (px >= 0 && py >= 0 && px < m_buffer.getWidth() && py < m_buffer.getHeight()) {
          PixelPosition pixelpos = {px, py};
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

void Renderer::drawCircle(const Vector& center, real_t radius, const Color& color)
{
  PixelPosition centerPx = worldToPixel(center);
  int_t radiusPx = static_cast<int_t>(std::round(radius * m_scale));

  for (int_t y = -radiusPx; y <= radiusPx; y++) {
    for (int_t x = -radiusPx; x <= radiusPx; x++) {
      if (x * x + y * y <= radiusPx * radiusPx) {
        int_t px = centerPx.x + x;
        int_t py = centerPx.y + y;
        if (px >= 0 && py >= 0 && px < m_buffer.getWidth() && py < m_buffer.getHeight()) {
          PixelPosition pixelpos = {px, py};
          m_buffer[pixelpos] = color;
        }
      }
    }
  }
}

void Renderer::drawPolygon(const std::vector<Vector>& vertices, uint_t widthPixels, const Color& color)
{
  if (vertices.size() < 3) {
    return;
  }
  for (size_t i = 0; i < vertices.size(); i++) {
    size_t next = (i + 1) % vertices.size();
    drawLine(vertices[i], vertices[next], widthPixels, color);
  }
}

void Renderer::drawArrow(const Vector& start, const Vector& end, uint_t widthPixels, const Color& color)
{
  drawLine(start, end, widthPixels, color);
  // TODO: Implement arrowhead drawing
}

void Renderer::drawGrid(uint_t widthPixels, const Color& color)
{
  real_t worldWidth = static_cast<real_t>(m_width) / m_scale;
  real_t worldHeight = static_cast<real_t>(m_height) / m_scale;

  int_t numLinesX = static_cast<int_t>(std::ceil(worldWidth / 2.0));
  int_t numLinesY = static_cast<int_t>(std::ceil(worldHeight / 2.0));

  // Draw vertical lines (parallel to Y-axis)
  for (int_t i = -numLinesX; i <= numLinesX; i++) {
    real_t x = static_cast<real_t>(i);
    Vector top = {x, worldHeight / 2.0};
    Vector bottom = {x, -worldHeight / 2.0};
    drawLine(top, bottom, widthPixels, color);
  }

  // Draw horizontal lines (parallel to X-axis)
  for (int_t j = -numLinesY; j <= numLinesY; j++) {
    real_t y = static_cast<real_t>(j);
    Vector left = {-worldWidth / 2.0, y};
    Vector right = {worldWidth / 2.0, y};
    drawLine(left, right, widthPixels, color);
  }
}

void Renderer::drawWorldFrame(uint_t widthPixels)
{
  Vector origin = {0.0, 0.0};
  Vector xAxis = {1.0, 0.0};
  Vector yAxis = {0.0, 1.0};
  drawLine(origin, xAxis, widthPixels, red);
  drawLine(origin, yAxis, widthPixels, blue);
}

void Renderer::display()
{
  bool shouldClose = false;
  while (!shouldClose) {
    // Update texture with framebuffer
    SDL_UpdateTexture(m_texture, nullptr, m_buffer.getRawData(), static_cast<int>(m_buffer.getWidth() * 4));
    SDL_SetRenderDrawColor(m_SDLRenderer, 0, 0, 0, 255);
    SDL_RenderClear(m_SDLRenderer);
    SDL_RenderCopy(m_SDLRenderer, m_texture, nullptr, nullptr);
    SDL_RenderPresent(m_SDLRenderer);

    // Handle events
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      if (event.type == SDL_QUIT || (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_ESCAPE)) {
        shouldClose = true;
      }
    }
    SDL_Delay(16);  // ~60 FPS
  }
}

#pragma endregion

}  // namespace sim