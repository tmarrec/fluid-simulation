#pragma once
#include <iostream>
#include <numeric>

#include "Renderer.h"
#include "Window.h"
#include "types.h"
#include "Components.h"
#include "Shader.h"
#include "ecs/Coordinator.h"
#include "BasicEntities.h"

#include "systems/MeshRenderer.h"
#include "systems/Fluids.h"

class Game
{
public:
	void run(WindowInfos windowInfos);

private:
	void mainLoop();
    void initECS();

	std::shared_ptr<Window> _window = std::make_shared<Window>();
	Renderer _renderer {};
    std::vector<Entity> _entities {};

    std::shared_ptr<MeshRenderer> _meshRendererSys = nullptr;
    std::shared_ptr<Fluids> _fluidsSys = nullptr;
};
