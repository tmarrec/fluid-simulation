#include "utils.h"
#include "config.h"
#include "Simulation.h"

int main()
{
    PRINT_TITLE();

    WindowInfos windowInfos;
    windowInfos.title = "fluid-simulation - Tristan Marrec";
    windowInfos.x = 800;
    windowInfos.y = 800;

    readConfig();

	Simulation sim;
	sim.run(windowInfos);
	return EXIT_SUCCESS;
}
