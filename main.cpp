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
  const real_t ANGLE = M_PI_2;
  Rotation<2> R2(ANGLE);
  Vector v2 = {1, 0};

  const int_t WIDTH = 1600;
  const int_t HEIGHT = 900;
  const real_t SCALE = 100.0;  // 100 pixels per unit

  Renderer rend(WIDTH, HEIGHT, SCALE);
  rend.clear(lightBlue);
  rend.drawGrid(1);
  rend.drawWorldFrame();

  // Example: draw a circle at world position (2, 3) with radius 0.5 units
  // rend.drawCircle({2.0, 3.0}, 0.5, green);

  rend.display();

  return 0;
}

int demo(void)
{
  const real_t ANGLE = M_PI_2;
  Vector v3 = {1, 0, 0};
  std::cout << Rotation<3>::Rz(ANGLE) * v3 << std::endl;
  return 0;
}