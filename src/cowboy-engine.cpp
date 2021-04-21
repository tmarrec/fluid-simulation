#include <iostream>

#include "types.h"
#include "utils.h"
#include "Game.h"

int main()
{
	PRINT_TITLE();

    WindowInfos windowInfos;
    windowInfos.title = "cowboy-engine";
    windowInfos.x = 1400;
    windowInfos.y = 800;

	Game game;
	game.run(windowInfos);
	return EXIT_SUCCESS;
}
