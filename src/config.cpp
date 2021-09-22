#include "config.h"

namespace Config
{
    std::uint16_t N = 64;
    std::uint16_t dim = 2;
    Solver solver = PCG;
    Advection advection = SEMI_LAGRANGIAN;
    double dt = 0.000004;
    bool exportFrames = false;
    bool renderFrames = true;
    std::uint16_t width = 800;
    std::uint16_t height = 800;
    std::uint64_t endFrame = 65536;
}  // namespace Config

// Read config.ini file and set variables of the global namespace Config
void readConfig()
{
    inipp::Ini<char> ini;
    std::ifstream is("config.ini");
    if (!is.is_open())
    {
        WARNING("No config.ini file! Will use default configuration values");
    }
    else
    {
        ini.parse(is);
        inipp::get_value(ini.sections["GRID"], "N",
                Config::N);
        inipp::get_value(ini.sections["GRID"], "dim",
                Config::dim);
        inipp::get_value(ini.sections["FLUID"], "dt",
                Config::dt);
        inipp::get_value(ini.sections["RENDER"], "exportFrames",
                Config::exportFrames);
        inipp::get_value(ini.sections["RENDER"], "renderFrames",
                Config::renderFrames);
        inipp::get_value(ini.sections["RENDER"], "width",
                Config::width);
        inipp::get_value(ini.sections["RENDER"], "height",
                Config::height);
        inipp::get_value(ini.sections["RENDER"], "endFrame",
                Config::endFrame);

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
}

