#include "config.h"

namespace Config
{
    std::uint16_t N = 256;
    std::uint16_t dim = 2;
    Solver solver = PCG;
    Advection advection = SEMI_LAGRANGIAN;
    double viscosity = 1.15;
    double diffusion = 1.00;
    double dt = 0.000004;
    bool exportFrames = false;
}

void readConfig()
{
    inipp::Ini<char> ini;
    std::ifstream is("config.ini");
    if (!is.is_open())
    {
        WARNING("No config file! Will use default configuration values");
    }
    else
    {
        ini.parse(is);
        inipp::get_value(ini.sections["GRID"], "N", Config::N);
        inipp::get_value(ini.sections["GRID"], "dim", Config::dim);
        inipp::get_value(ini.sections["FLUID"], "dt", Config::dt);
        inipp::get_value(ini.sections["FLUID"], "viscosity", Config::viscosity);
        inipp::get_value(ini.sections["FLUID"], "diffusion", Config::diffusion);
        inipp::get_value(ini.sections["RENDER"], "exportFrames", Config::exportFrames);

        if (!(Config::dim == 2 || Config::dim == 3))
        {
            ERROR("dim should be either 2 or 3");
        }

        std::string temp;

        inipp::get_value(ini.sections["SOLVER"], "solver", temp);
        if (temp == "CG")
            Config::solver = CG;
        else if (temp == "PCG")
            Config::solver = PCG;

        inipp::get_value(ini.sections["SOLVER"], "advection", temp);
        if (temp == "SEMI_LAGRANGIAN")
            Config::advection = SEMI_LAGRANGIAN;
        else if (temp == "MACCORMACK")
            Config::advection = MACCORMACK;
    }

    INFO("\033[1m=== CONFIGURATION ===\033[0m");
    INFO("\033[1m[GRID]\033[0m")
    INFO("N             = " << Config::N);
    INFO("dim           = " << Config::dim);
    INFO("\033[1m[SOLVER]\033[0m")
    INFO("solver        = " << Config::solver);
    INFO("advection     = " << Config::advection);
    INFO("\033[1m[FLUID]\033[0m")
    INFO("dt            = " << Config::dt);
    INFO("viscosity     = " << Config::viscosity);
    INFO("diffusion     = " << Config::diffusion);
    INFO("\033[1m[RENDER]\033[0m")
    INFO("exportFrames  = " << Config::exportFrames);
    INFO("\033[1m=====================\033[0m");
}

