#pragma once

#include "../ecs/Coordinator.h"
#include "../Components.h"
#include "../BasicEntities.h"
#include <cstdint>
#include <vector>
#include <unistd.h>
#include <ctime>

extern Coordinator gCoordinator;

class Fluids : public System
{
public:
    void init(std::shared_ptr<Renderer> renderer);
    void update();
    void reset(bool force = false);

#ifdef DEBUG_GUI
    void fluidDebugTool();
    void fluidSetupDebug();
#endif

private:
    void Vstep(Fluid3D& fluid);
    void Sstep(Fluid3D& fluid);

    void addSource(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& S) const;
    void diffuse(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& Xprev, float diff, std::uint8_t b) const;

    void advect(Fluid3D& fluid, std::vector<float>& D, std::vector<float>& Dprev, std::vector<float>& X, std::vector<float>& Y, std::vector<float>& Z, std::uint8_t b) const;
    void project(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& Y, std::vector<float>& Z, std::vector<float>& p, std::vector<float>& div) const;
    void GaussSeidelRelaxationLinSolve(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& Xprev, float a, float c, std::uint8_t b) const;

    void ConjugateGradientMethodLinSolve(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& Xprev, float a, float c, std::uint8_t b) const;

    void setBnd(Fluid3D& fluid, std::vector<float>& X, std::uint8_t b) const;

    void updateRender(Fluid3D& fluid);

    std::shared_ptr<Renderer> _renderer = nullptr;

#ifdef DEBUG_GUI
    float _debugViscosity;
    float _debugDiffusion;
    float _debugDt;
    float _debugAbsorption;
    float _debugLightIntensity[3];
    int _debugN;

    std::deque<float> _debugVstepTimes;
    std::deque<float> _debugSstepTimes;
    std::deque<float> _debugTextureTimes;

    std::deque<float> _debugVstepDiffuseTimes;
    std::deque<float> _debugVstepProjectTimes;
    std::deque<float> _debugVstepAdvectTimes;

    std::deque<float> _debugSstepDiffuseTimes;
    std::deque<float> _debugSstepAdvectTimes;
#endif

};

