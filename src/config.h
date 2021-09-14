#pragma once

#include "inipp.h"
#include "utils.h"
#include "types.h"

#include <fstream>

namespace Config
{
    extern std::uint16_t N;
    extern std::uint16_t dim;
    extern double viscosity;
    extern double diffusion;
    extern double dt;
    extern Solver solver;
    extern Advection advection;
    extern bool exportFrames;
}

void readConfig();
