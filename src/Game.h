#pragma once
#include <iostream>
#include <chrono>

#include "Renderer.h"
#include "Window.h"
#include "types.h"
#include "Components.h"
#include "Shader.h"
#include "ecs/Coordinator.h"
#include "BasicEntities.h"

#include "systems/Physics.h"
#include "systems/MeshRenderer.h"

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
    std::shared_ptr<MeshRenderer> _meshRendererSys;
};
