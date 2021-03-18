#pragma once
#include <iostream>

#include "Renderer.h"
#include "Window.h"

class Game
{
public:
	void run();

private:
	void renderInit();
	void windowInit();
	void mainLoop();

	std::shared_ptr<Window> _window = std::make_shared<Window>();
	Renderer _renderer;
};
