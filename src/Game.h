#pragma once
#include <iostream>

#include "Renderer.h"
#include "Window.h"
#include "ecs/System.h"
#include "systems/Physics.h"
#include "types.h"
#include "ecs/Coordinator.h"
#include "Components.h"


class Game
{
public:
	void run(WindowInfos windowInfos);

private:
	void mainLoop();

	std::shared_ptr<Window> _window = std::make_shared<Window>();
	Renderer _renderer;
    std::vector<Entity> _entities;
    std::shared_ptr<Physics> _physicsSys;
};
