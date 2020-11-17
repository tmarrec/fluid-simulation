#include "Game.h"
#include "Renderer.h"

void Game::run()
{
	windowInit();
	renderInit();
	mainLoop();
}

void Game::windowInit()
{
	_window->init();
}

void Game::renderInit()
{
	_renderer.init(_window);
}

void Game::mainLoop()
{
	while (!_window->windowShouldClose())
	{
		_window->pollEvents();
	}
}

