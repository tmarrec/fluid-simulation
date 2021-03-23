#include "Game.h"
#include <chrono>
#include <memory>

Coordinator gCoordinator;

void Game::run(WindowInfos windowInfos)
{
	_window->init(windowInfos);
	_renderer.init(_window);

    gCoordinator.Init();
    _entities.reserve(MAX_ENTITIES);

    gCoordinator.RegisterComponent<Transform>();
    gCoordinator.RegisterComponent<Mesh>();
    gCoordinator.RegisterComponent<Camera>();
    gCoordinator.RegisterComponent<Material>();

    _physicsSys = gCoordinator.RegisterSystem<Physics>();
    _meshRendererSys = gCoordinator.RegisterSystem<MeshRenderer>();
    Camera camera
    {
        .yaw = 0,
        .pitch = 0,
        .speed = 5,
        .FOV = 60
    };
    _meshRendererSys->init(std::make_shared<Renderer>(_renderer), camera);

    Signature signaturePhysics;
    signaturePhysics.set(gCoordinator.GetComponentType<Transform>());
    gCoordinator.SetSystemSignature<Physics>(signaturePhysics);

    Signature signatureMeshRenderer;
    signatureMeshRenderer.set(gCoordinator.GetComponentType<Transform>());
    signatureMeshRenderer.set(gCoordinator.GetComponentType<Mesh>());
    signatureMeshRenderer.set(gCoordinator.GetComponentType<Material>());
    gCoordinator.SetSystemSignature<MeshRenderer>(signatureMeshRenderer);

    auto entity = gCoordinator.CreateEntity();
    gCoordinator.AddComponent(entity, Transform
    {
        .position = glm::vec3(0, 0, 0),
        .rotation = glm::vec3(0, 0, 0),
        .scale = glm::vec3(1, 1, 1)
    });

    gCoordinator.AddComponent(entity, Mesh
    {
        .vertices =
        {
            -0.5f, -0.5f, 0.5f,
            0.5f, -0.5f, 0.5f,
            0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, 0.5f,

            0.5f, 0.5f, 0.5f,
            0.5f, 0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, 0.5f,

            -0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, 0.5f, -0.5f,
            -0.5f, 0.5f, -0.5f,

            -0.5f, -0.5f, -0.5f,
            -0.5f, -0.5f, 0.5f,
            -0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, -0.5f,
                            
            0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, 0.5f,
            -0.5f, 0.5f, -0.5f,
            0.5f, 0.5f, -0.5f,

            -0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, 0.5f,
            -0.5f, -0.5f, 0.5f,
        },
        .normals =
        {
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 1.0f,

            1.0f, 0.0f, 0.0f,
            1.0f, 0.0f, 0.0f,
            1.0f, 0.0f, 0.0f,
            1.0f, 0.0f, 0.0f,

            0.0f, 0.0f, -1.0f,
            0.0f, 0.0f, -1.0f,
            0.0f, 0.0f, -1.0f,
            0.0f, 0.0f, -1.0f,

            -1.0f, 0.0f, 0.0f,
            -1.0f, 0.0f, 0.0f,
            -1.0f, 0.0f, 0.0f,
            -1.0f, 0.0f, 0.0f,

            0.0f, 1.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 1.0f, 0.0f,

            0.0f, -1.0f, 0.0f,
            0.0f, -1.0f, 0.0f,
            0.0f, -1.0f, 0.0f,
            0.0f, -1.0f, 0.0f,
        },
        .indices =
        {
            0,  1,  2,  0,  2,  3,   // Front
            3,  6,  5,  4,  7,  6,   // Right
            5,  9, 11,  11, 9,  8,    // Back
            15, 12, 13, 13, 14, 15,  // Left
            17, 16, 19, 17, 19, 18,  // Upper
            20, 21, 23, 21, 22, 23 	 // Bottom
        },
    });

    Shader shaderProgram {};
    shaderProgram.setVert("shaders/vert.vert");
    shaderProgram.setFrag("shaders/frag.frag");

    gCoordinator.AddComponent(entity, Material
    {
        .shader = shaderProgram
    });

	mainLoop();
}

void Game::mainLoop()
{
    float dt = 0.0f;
	while (!_window->windowShouldClose())
	{
        auto startTime = std::chrono::high_resolution_clock::now();

        _physicsSys->update(dt);
        _renderer.prePass();
        _meshRendererSys->update();

        _window->swapBuffers();
		_window->pollEvents();

        auto stopTime = std::chrono::high_resolution_clock::now();
        dt = std::chrono::duration<float, std::chrono::seconds::period>(stopTime - startTime).count();
        //std::cout << 1/dt << std::endl;
	}
}

