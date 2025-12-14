/*****************************************************************************************
 * @file: main.cpp
 *
 * @brief: 3D Rendering Demo
 *
 * @details: Demonstrates 3D rendering with rigid bodies including rotating shapes.
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
#include "Material.hpp"
#include "Mesh.hpp"
#include "Renderer.hpp"
#include "RigidBody.hpp"
#include "Rotation.hpp"
#include "Transform.hpp"
#include "Vector.hpp"
#include "types.hpp"

using namespace sim;

int main(void)
{
    // Create camera looking at the scene
    // Camera positioned to look directly at objects
    Transform cam_transform(Rotation(), Vector({0, 0, 2}));
    Camera camera(1600, 900, cam_transform, deg2rad(60), 0.1, 100);

    // Create renderer with camera
    Renderer renderer(camera, "3D Scene Demo - Rigid Bodies");

    // Create materials
    Material box_material(1000.0, Color::Cyan());
    Material sphere_material(800.0, Color::Magenta());

    // Position vectors
    Vector box_position = {1, 0, 0};
    Vector sphere_position = {-1, 0, 0};

    // Create rigid bodies with meshes
    RigidBody box(Mesh::Cube(0.25, box_material.color), box_material, Transform(box_position));

    RigidBody sphere(Mesh::Sphere(0.25, sphere_material.color), sphere_material, Transform(sphere_position));

    real_t angle = 0.0;

    // Main rendering loop
    while (!renderer.shouldClose())
    {
        renderer.pollEvents();

        // Update rotation
        angle += 0.02;

        // Update transforms of rotating bodies
        box.setTransform(Transform(Rotation(angle * 0.3, angle, angle * 0.7), box_position));

        sphere.setTransform(Transform(Rotation(angle * 0.5, angle * 1.3, 0.0), sphere_position));

        // Clear framebuffer
        renderer.clear(Color(0.1, 0.1, 0.15, 1.0));  // Dark blue background

        // Draw world coordinate axes
        renderer.renderWorldFrame(1.0);

        // Render all rigid bodies
        renderer.renderBody(box);
        renderer.renderBody(sphere);

        // Present to screen
        renderer.present();

        // Frame rate control (~60 FPS)
        SDL_Delay(16);
    }

    return 0;
}
