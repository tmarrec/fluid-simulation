#include "Game.h"
#include "BasicEntities.h"

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
        .speed = 0.1f,
        .FOV = 90,
        .transform = Transform
            {
                .position = {0, 0, -2},
                .rotation = {0, 0, 0},
                .scale = {1, 1, 1}
            }
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

    BasicEntities::addCube();
    
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

