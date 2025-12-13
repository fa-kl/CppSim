/*****************************************************************************************
 * @file: Mesh.cpp
 *
 * @brief: This file implements a mesh class for rendering.
 *
 * @author: fakl
 * @date: December 2025
 *
 ****************************************************************************************/

#include "Mesh.hpp"

namespace sim
{

Mesh operator*(const Transform& transform, const Mesh& mesh)
{
    Mesh result = mesh;
    for (Triangle& triangle : result)
    {
        for (Vertex& vert : triangle)
        {
            vert.point = transform * vert.point;
        }
    }
    return result;
}

Mesh cullMesh(const Mesh& input_mesh, const Camera& camera)
{
    Mesh result;
    result.reserve(input_mesh.size());
    for (const Triangle& triangle : input_mesh)
    {
        bool any_in_view = false;
        for (const Vertex& vertex : triangle)
        {
            if (camera.isInView(vertex.point))
            {
                any_in_view = true;
                break;
            }
        }
        if (any_in_view)
        {
            result.push_back(triangle);
        }
    }
    return result;
}

}  // namespace sim