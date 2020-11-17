#include <iostream>

#define GLFW_INCLUDE_VULKAN
#include "utils.h"
#include "Game.h"

int main()
{
	PRINT_TITLE();
	Game game;
	game.run();
	return EXIT_SUCCESS;
}