#include "Game.h"
#include "Renderer.h"

void Game::run()
{
	renderInit();
	windowInit();
	mainLoop();
}

void Game::renderInit()
{
	renderer.init();
}

void Game::windowInit()
{
	window.init();
}

void Game::mainLoop()
{
	while (!window.windowShouldClose())
	{
		window.pollEvents();
	}
}

