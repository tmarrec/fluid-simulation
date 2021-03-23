#include "Game.h"
#include <chrono>

Coordinator gCoordinator;

void Game::run(WindowInfos windowInfos)
{
	_window->init(windowInfos);
	_renderer.init(_window);

    gCoordinator.Init();
    _entities.reserve(MAX_ENTITIES);

    gCoordinator.RegisterComponent<Transform>();
    gCoordinator.RegisterComponent<Mesh>();

    _physicsSys = gCoordinator.RegisterSystem<Physics>();
    _meshRendererSys = gCoordinator.RegisterSystem<MeshRenderer>();
    _meshRendererSys->init(std::make_shared<Renderer>(_renderer));

    Signature signaturePhysics;
    signaturePhysics.set(gCoordinator.GetComponentType<Transform>());
    gCoordinator.SetSystemSignature<Physics>(signaturePhysics);

    Signature signatureMeshRenderer;
    signatureMeshRenderer.set(gCoordinator.GetComponentType<Transform>());
    signatureMeshRenderer.set(gCoordinator.GetComponentType<Mesh>());
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
             0.5f,  0.5f, 0.0f,  // top right
             0.5f, -0.5f, 0.0f,  // bottom right
            -0.5f, -0.5f, 0.0f,  // bottom left
            -0.5f,  0.5f, 0.0f   // top left 
        },
        .normals =
        {
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
        },
        .indices =
        {
            0, 1, 3,
            1, 2, 3
        },
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

