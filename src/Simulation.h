#pragma once

#include <iostream>
#include <numeric>
#include <vector>

#include "./types.h"
#include "./Renderer.h"
#include "./Shader.h"
#include "./Window.h"
#include "./Fluids.h"
#include "./Input.h"
#include "./MarchingCube.h"

class Simulation
{
 public:
    void initRendering();
    void run();

 private:
    void initSimulationRendering();
    void updateMeshVec();
    void updateMeshGrid();
    void updateMeshGridBorder();
    void setCameraDir();
    void handleInputs();
    void stepFluid(const std::uint64_t it);
    void renderFrame();
    void printStatus(const std::uint64_t it, const float dt) const;

    Window _window = {};
    Renderer _renderer = {};
    MarchingCube marchingCube;

    Camera _camera = {};
    Fluids _fluid;
    FluidRenderer _fluidRenderer =
    {
        .transform = Transform
        {
            .position = {0, 0, 0},
            .rotation = {0, 0, 0},
            .scale = {1, 1, 1}
        },
        .mesh = Mesh
        {
            .vertices =
            {
                -0.5f, 0.0f, -0.5f,
                0.5f, 0.0f, -0.5f,
                0.5f, 0.0f, 0.5f,
                -0.5f, 0.0f, 0.5f,
            },
            .normals =
            {
                0.0f, 0.0f, 1.0f,
                0.0f, 0.0f, 1.0f,
                0.0f, 0.0f, 1.0f,
                0.0f, 0.0f, 1.0f,
            },
            .indices =
            {
                0,  1,  2,  0,  2,  3,
            },
        },
        .meshVec = Mesh
        {
            .vertices = {},
            .normals = {},
            .indices = {},
            .renderMode = LINES,
        },
        .meshGrid = Mesh
        {
            .vertices = {},
            .normals = {},
            .indices = {},
            .renderMode = LINES,
        },
        .meshGridBorder = Mesh
        {
            .vertices = {},
            .normals = {},
            .indices = {},
            .renderMode = LINES,
        },
        .material = {},
        .materialVec = {},
        .materialGrid = {},
        .materialGridBorder = {},
    };
};
