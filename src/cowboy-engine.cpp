#include "utils.h"
#include "Simulation.h"

int main()
{
    PRINT_TITLE();

    WindowInfos windowInfos;
    windowInfos.title = "cowboy-engine";
    windowInfos.x = 1000;
    windowInfos.y = 1000;

	Simulation sim;
	sim.run(windowInfos);
	return EXIT_SUCCESS;
}
