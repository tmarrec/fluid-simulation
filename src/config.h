#pragma once

#include <fstream>
#include <string>

#include "./inipp.h"
#include "./utils.h"
#include "./types.h"

namespace Config
{
    extern std::uint16_t N;
    extern std::uint16_t dim;
    extern double dt;
    extern Solver solver;
    extern Advection advection;
    extern bool exportFrames;
    extern bool renderFrames;
    extern std::uint16_t width;
    extern std::uint16_t height;
    extern std::uint64_t endFrame;
}  // namespace Config

void readConfig();
