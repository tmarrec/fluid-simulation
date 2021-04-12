#include "Fluids.h"
#include <cstdint>
#include <unistd.h>

void Fluids::update()
{
    for (auto const& entity : mEntities)
    {
        auto& fluid = gCoordinator.GetComponent<Fluid2D>(entity);
        auto& transform = gCoordinator.GetComponent<Transform>(entity);
        Vstep(fluid);
        Sstep(fluid);
        updateRender(fluid, transform);
    }
}

void Fluids::updateRender(Fluid2D& fluid, Transform& transform)
{
    // Testing purpose
    for (std::uint32_t i = 0; i < (fluid.N+2)*(fluid.N+2); ++i)
    {
        auto& transform = gCoordinator.GetComponent<Transform>(fluid.entities[i]);
        transform.rotation = fluid.U[i];
    }
}

void Fluids::Vstep(Fluid2D& fluid)
{
    addSource(fluid);
    fluid.Uprev = fluid.U;
    diffuse(fluid, 1.0f, fluid.U, fluid.Uprev);
    fluid.Uprev = fluid.U;
    //advect();
    usleep(20000);
}

void Fluids::Sstep(Fluid2D& fluid)
{
    /*
    addSource(fluid);
    fluid.Uprev = fluid.U;
    diffuse(fluid, 1.0f, fluid.U, fluid.Uprev);
    fluid.Uprev = fluid.U;
    //advect();
    usleep(20000);
    */
}

void Fluids::addSource(Fluid2D& fluid)
{
    // Testing
    fluid.U[44] = glm::vec3{0,4500.0f,0};
}

void Fluids::diffuse(Fluid2D& fluid, float visc, std::vector<glm::vec3>& X, std::vector<glm::vec3>& X0)
{
    const int kMax = 1;
    float a = fluid.dt * visc * fluid.N * fluid.N;

    for (std::uint8_t k = 0; k < kMax; ++k)
    {
        for (std::uint32_t x = 1; x <= fluid.N; ++x)
        {
            for (std::uint32_t y = 1; y <= fluid.N; ++y)
            {
                X[fluid.IX(x,y)] = 
                    (X0[fluid.IX(x,y)] + 
                        a*(X[fluid.IX(x-1,y)]+X[fluid.IX(x+1,y)]+X[fluid.IX(x,y-1)]+X[fluid.IX(x,y+1)]))
                    /(1+4*a);
            }
        }
        setBnd(fluid, X);
    }
}

void Fluids::advect(Fluid2D& fluid, std::vector<glm::vec3>& U, std::vector<glm::vec3>& Uprev)
{
}

void Fluids::setBnd(Fluid2D& fluid, std::vector<glm::vec3>& U)
{
    for (std::uint32_t i = 0; i <= fluid.N+1; ++i)
    {
        U[fluid.IX(0,i)] = glm::vec3{0,0,0};
        U[fluid.IX(i,0)] = glm::vec3{0,0,0};
        U[fluid.IX(fluid.N+1,i)] = glm::vec3{0,0,0};
        U[fluid.IX(i,fluid.N+1)] = glm::vec3{0,0,0};
    }
}
