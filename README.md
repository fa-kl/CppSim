# CppSim
A physics engine in C++

This is a project with the aim of gaining knowledge about simulators or physics engines with the goal of being able to build and control some simple mechanisms like a (double) pendulum on a cart.

The project is still in its early stages... stay tuned!

### Build Instructions (Ubuntu)
Install the required packages and build the project:
```bash
sudo apt install -y build-essential cmake libsdl2-dev
mkdir build && cd build && cmake .. && make
```

### Running the 3D Renderer Demo
After building, run the demo to see a rotating 3D box with rasterized triangles:
```bash
cd build
./Main
```

The demo features:
- Real-time 3D rendering with perspective projection
- Rotating box with colored fill and wireframe edges
- View frustum culling
- Press ESC to quit
