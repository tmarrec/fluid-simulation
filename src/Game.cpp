#include "Game.h"
#include "Components.h"
#include "ecs/Coordinator.h"
#include "systems/Physics.h"
#include "types.h"

Coordinator gCoordinator;

void Game::run(WindowInfos windowInfos)
{
	_window->init(windowInfos);
	_renderer.init(_window);

    gCoordinator.Init();
    _entities.reserve(MAX_ENTITIES);

    gCoordinator.RegisterComponent<Transform>();

    _physicsSys = gCoordinator.RegisterSystem<Physics>();

    Signature signature;
    signature.set(gCoordinator.GetComponentType<Transform>());
    gCoordinator.SetSystemSignature<Physics>(signature);

    auto entity = gCoordinator.CreateEntity();
    gCoordinator.AddComponent(entity, Transform
    {
        .position = glm::vec3(0, 0, 0),
        .rotation = glm::vec3(0, 0, 0),
        .scale = glm::vec3(1, 1, 1)
    });

	mainLoop();
}

void Game::mainLoop()
{
	while (!_window->windowShouldClose())
	{
        _physicsSys->Update(0.1f);

        _renderer.pass();
        _window->swapBuffers();
		_window->pollEvents();
	}
}

