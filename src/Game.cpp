#include "Game.h"

void Game::run(WindowInfos windowInfos)
{
	windowInit(windowInfos);
	renderInit();
	mainLoop();
}

void Game::windowInit(WindowInfos windowInfos)
{
	_window->init(windowInfos);
}

void Game::renderInit()
{
	_renderer.init(_window);
}

void Game::mainLoop()
{
	while (!_window->windowShouldClose())
	{
        _renderer.pass();
        _window->swapBuffers();
		_window->pollEvents();
	}
}

