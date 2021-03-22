#pragma once
#include <iostream>

#include "Renderer.h"
#include "Window.h"
#include "types.h"

class Game
{
public:
	void run(WindowInfos windowInfos);

private:
	void windowInit(WindowInfos windowInfos);
	void renderInit();
	void mainLoop();

	std::shared_ptr<Window> _window = std::make_shared<Window>();
	Renderer _renderer;
};
