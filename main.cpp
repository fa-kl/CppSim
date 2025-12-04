/*****************************************************************************************
 * @file: main.cpp
 *
 * @brief: This is the simulator's file.
 *
 * @details: This file may be used to demonstrate the simulator's functionality.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include <cmath>
#include <iostream>
#include <vector>

#include "Matrix.hpp"
#include "Renderer.hpp"
#include "Rotation.hpp"
#include "Shapes.hpp"
#include "Transform.hpp"
#include "Vector.hpp"
#include "types.hpp"

/*****************************************************************************************
 * Main Function
 ****************************************************************************************/

using namespace sim;

int main(void)
{
  try {
    const uint_t WIDTH = 800;
    const uint_t HEIGHT = 600;

    // Create window and framebuffer
    Window window(WIDTH, HEIGHT);
    FrameBuffer frameBuffer(WIDTH, HEIGHT);
    Renderer renderer(frameBuffer);

    // Clear with dark background
    frameBuffer.clear(black);

    // Draw grid
    /*
    const uint_t d = 100;
    for (uint_t x = d; x < WIDTH; x += d) {
      renderer.drawLine({x, 0}, {x, HEIGHT - 1}, 1, green);
    }
    for (uint_t y = d; y < HEIGHT; y += d) {
      renderer.drawLine({0, y}, {WIDTH - 1, y}, 1, blue);
    }
    */

    // Main rendering loop
    while (!window.shouldClose()) {
      window.display(frameBuffer);
      SDL_Delay(16);  // ~60 FPS
    }

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  std::cout << "Feel free to try out some stuff here ..." << std::endl;
  return 0;
}
