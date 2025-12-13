/*****************************************************************************************
 * @file: main.cpp
 *
 * @brief: 3D Rendering Demo
 *
 * @details: Demonstrates 3D rendering with a rotating box, sphere, floor, and world axes.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include <cmath>
#include <iostream>

#ifdef __linux__
#include <SDL2/SDL.h>
#else
#include <SDL.h>
#endif

#include "Camera.hpp"
#include "Color.hpp"
#include "Renderer.hpp"
#include "Rotation.hpp"
#include "Shapes.hpp"
#include "Transform.hpp"
#include "Vector.hpp"
#include "types.hpp"

using namespace sim;

int main(void)
{
    // Create renderer
    Renderer renderer(1600, 900, "3D Scene Demo");

    // Create camera looking at the scene
    Rotation cam_rot(deg2rad(-20), 0.0, 0.0);  // Slight tilt downward
    Transform cam_transform(cam_rot, Vector({3, 3, 10}));
    Camera camera(cam_transform, 800.0 / 600.0, deg2rad(60), 0.1, 100);

    // Create shapes
    Box box(1.0, CYAN);
    Sphere sphere(0.8, MAGENTA);
    Box floor(10.0, 0.2, 10.0, GRAY);

    real_t angle = 0.0;

    // Main rendering loop
    while (!renderer.shouldClose())
    {
        renderer.pollEvents();

        // Update rotation
        angle += 0.02;

        // Clear framebuffer
        renderer.clear(Color(0.1, 0.1, 0.15, 1.0));  // Dark blue background

        // Draw world coordinate axes
        renderer.drawAxes(camera, 2.0);

        // Render floor
        Rotation floor_rot(0.0, 0.0, 0.0);
        Transform floor_transform(floor_rot, Vector({0, -2, 0}));
        Mesh floor_mesh = floor_transform * floor.mesh();
        Mesh culled_floor = cullMesh(floor_mesh, camera);
        renderer.renderMesh(culled_floor, camera, DARK_GRAY);

        // Render rotating box
        Rotation box_rot(angle * 0.3, angle, angle * 0.7);
        Transform box_transform(box_rot, Vector({-1.5, 0, 0}));
        Mesh box_mesh = box_transform * box.mesh();
        Mesh culled_box = cullMesh(box_mesh, camera);
        renderer.renderMesh(culled_box, camera, WHITE);

        // Render rotating sphere
        Rotation sphere_rot(angle * 0.5, angle * 1.3, 0.0);
        Transform sphere_transform(sphere_rot, Vector({1.5, 0, 0}));
        Mesh sphere_mesh = sphere_transform * sphere.mesh();
        Mesh culled_sphere = cullMesh(sphere_mesh, camera);
        renderer.renderMesh(culled_sphere, camera, YELLOW);

        // Present to screen
        renderer.present();

        // Frame rate control (~60 FPS)
        SDL_Delay(16);
    }

    return 0;
}
