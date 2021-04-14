#pragma once

#include "../ecs/Coordinator.h"
#include "../Components.h"
#include "../BasicEntities.h"
#include <cstdint>

extern Coordinator gCoordinator;

struct Bound
{
    std::uint8_t dim;
    std::uint8_t b;
};

class Fluids : public System
{
public:
    void init(std::shared_ptr<Renderer> renderer);
    void update();

private:
    void Vstep(Fluid2D& fluid);
    void Sstep(Fluid2D& fluid);

    void addSource(Fluid2D& fluid, std::vector<float>& X, std::vector<float>& S);
    void diffuse(Fluid2D& fluid, std::vector<float>& X, std::vector<float>& Xprev, std::uint8_t b);
    void advect(Fluid2D& fluid, std::vector<float>& D, std::vector<float>& Dprev, std::vector<float>& X, std::vector<float>& Y, std::uint8_t b);
    void project(Fluid2D& fluid, std::vector<float>& X, std::vector<float>& Y, std::vector<float>& p, std::vector<float>& div);
    void setBnd(Fluid2D& fluid, std::vector<float>& X, std::uint8_t b);

    void updateRender(Fluid2D& fluid);

    std::shared_ptr<Renderer> _renderer = nullptr;
};

