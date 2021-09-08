#pragma once

#include "../inipp/inipp/inipp.h"
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
}

void readConfig();
