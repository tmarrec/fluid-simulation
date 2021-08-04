#include "utils.h"
#include "Simulation.h"

int main()
{
    PRINT_TITLE();

    WindowInfos windowInfos;
    windowInfos.title = "fluid-simulation - Tristan Marrec";
    windowInfos.x = 800;
    windowInfos.y = 800;

	Simulation sim;
	sim.run(windowInfos);
	return EXIT_SUCCESS;
}
