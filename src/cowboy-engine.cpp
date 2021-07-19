#include <iostream>

#include "types.h"
#include "utils.h"
#include "Game.h"

int main()
{
	//PRINT_TITLE();

    WindowInfos windowInfos;
    windowInfos.title = "cowboy-engine";
    windowInfos.x = 1000;
    windowInfos.y = 1000;

	Game game;
	game.run(windowInfos);
	return EXIT_SUCCESS;
}
