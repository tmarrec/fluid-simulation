#include "Game.h"
#include "BasicEntities.h"
#include <cstdint>

Coordinator gCoordinator;

void Game::initECS()
{
    gCoordinator.Init();
    _entities.reserve(MAX_ENTITIES);

    gCoordinator.RegisterComponent<Transform>();
    gCoordinator.RegisterComponent<Mesh>();
    gCoordinator.RegisterComponent<Camera>();
    gCoordinator.RegisterComponent<Material>();
    gCoordinator.RegisterComponent<Fluid3D>();

    _meshRendererSys = gCoordinator.RegisterSystem<MeshRenderer>();
    Camera camera
    {
        //.yaw = 225,
        //.pitch = -29,
        .yaw = -90,
        .pitch = -90,
        .speed = 0.1f,
        .FOV = 90,
        .transform = Transform
            {
                //.position = {9.0, 6.26, 9.0},
                .position = {0.0, 6.5, 0.0},
                .rotation = {0, 0, 0},
                .scale = {1, 1, 1}
            }
    };
    _meshRendererSys->init(_renderer, camera);
    _fluidsSys = gCoordinator.RegisterSystem<Fluids>();
    _fluidsSys->init(_renderer);

    Signature signatureMeshRenderer;
    signatureMeshRenderer.set(gCoordinator.GetComponentType<Transform>());
    signatureMeshRenderer.set(gCoordinator.GetComponentType<Mesh>());
    signatureMeshRenderer.set(gCoordinator.GetComponentType<Material>());
    gCoordinator.SetSystemSignature<MeshRenderer>(signatureMeshRenderer);

    Signature signatureFluids;
    signatureFluids.set(gCoordinator.GetComponentType<Transform>());
    signatureFluids.set(gCoordinator.GetComponentType<Fluid3D>());
    gCoordinator.SetSystemSignature<Fluids>(signatureFluids);
}

void Game::run(WindowInfos windowInfos)
{
	_window->init(windowInfos);
	_renderer->init(_window);
    initECS();

    BasicEntities::initBasicEntities(_renderer);

    BasicEntities::addFluid2D(glm::vec3{0,0,0}, glm::vec3{0,0,0}, glm::vec3{10,10,10});

    //BasicEntities::addLineCube(glm::vec3{5,0,0}, glm::vec3{0,0,0}, glm::vec3{5,5,5});

	mainLoop();
}

void Game::mainLoop()
{
    [[maybe_unused]] float dt = 0.0f;
    std::uint32_t iterations = 0;

	while (!_window->windowShouldClose())
	{
        if (iterations == 64)
            exit(0);
        //std::cout << "== Iteration " << iterations << " ==" << std::endl;
            
        auto startTime = std::chrono::high_resolution_clock::now();

        _fluidsSys->update(iterations);

        _renderer->prePass();
        _meshRendererSys->update();
        _renderer->endPass();

        _renderer->writeImg(iterations);
        _window->swapBuffers();

		_window->pollEvents();
        auto stopTime = std::chrono::high_resolution_clock::now();
        dt = std::chrono::duration<float, std::chrono::seconds::period>(stopTime - startTime).count();
        //std::cout << "Done in " << dt << " sec" << std::endl << std::endl;

        iterations++;
	}
}

