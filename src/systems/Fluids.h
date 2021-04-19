#pragma once

#include "../ecs/Coordinator.h"
#include "../Components.h"
#include "../BasicEntities.h"
#include <cstdint>
#include <vector>

extern Coordinator gCoordinator;

class Fluids : public System
{
public:
    void init(std::shared_ptr<Renderer> renderer);
    void update();

private:
    void Vstep(Fluid3D& fluid);
    void Sstep(Fluid3D& fluid);

    void addSource(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& S) const;
    void diffuse(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& Xprev, std::uint8_t b) const;
    void advect(Fluid3D& fluid, std::vector<float>& D, std::vector<float>& Dprev, std::vector<float>& X, std::vector<float>& Y, std::vector<float>& Z, std::uint8_t b) const;
    void project(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& Y, std::vector<float>& Z, std::vector<float>& p, std::vector<float>& div) const;
    void linSolve(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& Xprev, float a, float c, std::uint8_t b) const;
    void setBnd(Fluid3D& fluid, std::vector<float>& X, std::uint8_t b) const;
    void swap(std::vector<float>& X, std::vector<float>& Y) const;

    void updateRender(Fluid3D& fluid);

    std::shared_ptr<Renderer> _renderer = nullptr;
};

