#include "utils.h"
#include "config.h"
#include "Simulation.h"

int main()
{
    PRINT_TITLE();

    readConfig();

    Simulation sim;
    if (Config::renderFrames)
    {
        sim.initRendering();
    }
    sim.run();
    return EXIT_SUCCESS;
}
