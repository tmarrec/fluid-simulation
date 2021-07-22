#pragma once

#include "types.h"
#include "Renderer.h"
#include "Shader.h"
#include "Window.h"
#include "Fluids.h"

#include <iostream>
#include <numeric>

class Simulation
{
public:
	void run(WindowInfos windowInfos);

private:
	void mainLoop();
    void initSimulation();
    void updateMeshVec();

	Window _window = {};
    Renderer _renderer = {};

    Camera _camera = 
    {
        .yaw = -90,
        .pitch = -90,
        .speed = 0.1f,
        .FOV = 60,
        .transform = Transform
            {
                .position = {0.0, 1, 0.0},
                .rotation = {0, 0, 0},
                .scale = {1, 1, 1}
            }
    };
    Fluids _fluid = {};
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
            .is2D = true,
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
        .material = {},
        .materialVec = {},
    };
};
