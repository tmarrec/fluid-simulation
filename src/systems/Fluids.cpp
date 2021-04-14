#include "Fluids.h"
#include <cstdint>

void Fluids::init(std::shared_ptr<Renderer> renderer)
{
    _renderer = renderer;
}

void Fluids::update()
{
    for (auto const& entity : mEntities)
    {
        auto& fluid = gCoordinator.GetComponent<Fluid2D>(entity);

        fluid.substanceField[fluid.N*(fluid.N/2)+10] = 500;
        fluid.velocityFieldX[fluid.N*(fluid.N/2)+9] = 500;

        fluid.substanceField[fluid.N*(fluid.N/2)+(fluid.N-10)] = 500;
        fluid.velocityFieldX[fluid.N*(fluid.N/2)+(fluid.N-9)] = -500;

        Vstep(fluid);
        Sstep(fluid);
        updateRender(fluid);
    }
}

void Fluids::updateRender(Fluid2D& fluid)
{
    std::vector<std::uint8_t> texture((fluid.N+2)*(fluid.N+2)*3, 0);
    for (std::uint32_t i = 0; i < (fluid.N+2)*(fluid.N+2); ++i)
    {
        std::uint8_t density = std::clamp(static_cast<int>(fluid.substanceField[i]), 0, 255);
        texture[i*3] = density;
        texture[i*3+1] = density;
        texture[i*3+2] = density;
    }
    _renderer->updateTexture(texture);
}

void Fluids::Vstep(Fluid2D& fluid)
{
    addSource(fluid, fluid.velocityFieldX, fluid.velocityFieldPrevX);
    addSource(fluid, fluid.velocityFieldY, fluid.velocityFieldPrevY);

    fluid.velocityFieldPrevX = fluid.velocityFieldX;
    diffuse(fluid, fluid.velocityFieldX, fluid.velocityFieldPrevX, 1);
    fluid.velocityFieldPrevY = fluid.velocityFieldY;
    diffuse(fluid, fluid.velocityFieldY, fluid.velocityFieldPrevY, 2);

    project(fluid, fluid.velocityFieldX, fluid.velocityFieldY, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY);

    fluid.velocityFieldPrevX = fluid.velocityFieldX;
    fluid.velocityFieldPrevY = fluid.velocityFieldY;

    advect(fluid, fluid.velocityFieldX, fluid.velocityFieldPrevX, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, 1);
    advect(fluid, fluid.velocityFieldY, fluid.velocityFieldPrevY, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, 2);

    project(fluid, fluid.velocityFieldX, fluid.velocityFieldY, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY);
}

void Fluids::Sstep(Fluid2D& fluid)
{
    addSource(fluid, fluid.substanceField, fluid.substanceFieldPrev);

    fluid.substanceFieldPrev = fluid.substanceField;

    diffuse(fluid, fluid.substanceField, fluid.substanceFieldPrev, 0);

    fluid.substanceFieldPrev = fluid.substanceField;

    advect(fluid, fluid.substanceField, fluid.substanceFieldPrev, fluid.velocityFieldX, fluid.velocityFieldY, 0);
}

void Fluids::addSource(Fluid2D& fluid, std::vector<float>& X, std::vector<float>& S)
{
    for (std::uint32_t i = 0; i < (fluid.N+2)*(fluid.N+2); ++i)
    {
        X[i] += fluid.dt * S[i];
    }
}

void Fluids::diffuse(Fluid2D& fluid, std::vector<float>& X, std::vector<float>& Xprev, std::uint8_t b)
{
    const float a = fluid.dt * fluid.viscosity * fluid.N * fluid.N;

    for (std::uint8_t k = 0; k < 20; ++k)
    {
        for (std::uint32_t x = 1; x <= fluid.N; ++x)
        {
            for (std::uint32_t y = 1; y <= fluid.N; ++y)
            {
                X[fluid.IX(x,y)] = 
                    (Xprev[fluid.IX(x,y)] + 
                        a*(X[fluid.IX(x-1,y)]+X[fluid.IX(x+1,y)]+X[fluid.IX(x,y-1)]+X[fluid.IX(x,y+1)]))
                    /(1+4*a);
            }
        }
        setBnd(fluid, X, b);
    }
}

void Fluids::advect(Fluid2D& fluid, std::vector<float>& D, std::vector<float>& Dprev, std::vector<float>& X, std::vector<float>& Y, std::uint8_t b)
{
    float dt0 = fluid.dt * fluid.N;
    for (std::uint32_t i = 1; i <= fluid.N; ++i)
    {
        for (std::uint32_t j = 1; j <= fluid.N; ++j)
        {
            float x = i-dt0*X[fluid.IX(i,j)];
            float y = j-dt0*Y[fluid.IX(i,j)];

            x = std::clamp(x, 0.5f, fluid.N + 0.5f);
            std::uint32_t i0 = static_cast<int>(x);
            std::uint32_t i1 = i0 + 1;
            y = std::clamp(y, 0.5f, fluid.N + 0.5f);
            std::uint32_t j0 = static_cast<int>(y);
            std::uint32_t j1 = j0 + 1;

            float s1 = x    - i0;
            float s0 = 1.0f - s1;
            float t1 = y    - j0;
            float t0 = 1.0f - t1;

            D[fluid.IX(i,j)] = s0*(t0*Dprev[fluid.IX(i0,j0)]+t1*Dprev[fluid.IX(i0,j1)])+s1*(t0*Dprev[fluid.IX(i1,j0)]+t1*Dprev[fluid.IX(i1,j1)]);
        }
    }
    setBnd(fluid, D, b);
}

void Fluids::project(Fluid2D& fluid, std::vector<float>& X, std::vector<float>& Y, std::vector<float>& p, std::vector<float>& div)
{
    float h = 1.0f/fluid.N;
    for (std::uint32_t i = 1; i <= fluid.N; ++i)
    {
        for (std::uint32_t j = 1; j <= fluid.N; ++j)
        {
            div[fluid.IX(i,j)] = -0.5f*h*(X[fluid.IX(i+1,j)]-X[fluid.IX(i-1,j)]+Y[fluid.IX(i,j+1)]-Y[fluid.IX(i,j-1)]);
            p[fluid.IX(i,j)] = 0;
        }
    }
    setBnd(fluid, p, 0);
    setBnd(fluid, div, 0);

    for (std::uint32_t k = 0; k < 20; ++k)
    {
        for (std::uint32_t i = 1; i <= fluid.N; ++i)
        {
            for (std::uint32_t j = 1; j <= fluid.N; ++j)
            {
                p[fluid.IX(i,j)] = (div[fluid.IX(i,j)]+p[fluid.IX(i-1,j)]+p[fluid.IX(i+1,j)]+p[fluid.IX(i,j-1)]+p[fluid.IX(i,j+1)])/4;
            }
        }
        setBnd(fluid, p, 0);
    }
    
    for (std::uint32_t i = 1; i <= fluid.N; ++i)
    {
        for (std::uint32_t j = 1; j <= fluid.N; ++j)
        {
            X[fluid.IX(i,j)] -= 0.5f*(p[fluid.IX(i+1,j)]-p[fluid.IX(i-1,j)])/h;
            Y[fluid.IX(i,j)] -= 0.5f*(p[fluid.IX(i,j+1)]-p[fluid.IX(i,j-1)])/h;
        }
    }
    setBnd(fluid, X, 1);
    setBnd(fluid, Y, 2);
}

void Fluids::setBnd(Fluid2D& fluid, std::vector<float>& X, std::uint8_t b)
{
    for (std::uint32_t i = 1; i <= fluid.N; ++i)
    {
        X[fluid.IX(0,i)]         = b == 1 ? -X[fluid.IX(1,i)]       : X[fluid.IX(1,i)];
        X[fluid.IX(fluid.N+1,i)] = b == 1 ? -X[fluid.IX(fluid.N,i)] : X[fluid.IX(fluid.N,i)];
        X[fluid.IX(i,0)]         = b == 2 ? -X[fluid.IX(i,1)]       : X[fluid.IX(i,1)];
        X[fluid.IX(i,fluid.N+1)] = b == 2 ? -X[fluid.IX(i,fluid.N)] : X[fluid.IX(i,fluid.N)];
    }
    X[fluid.IX(0,0)]                    = 0.5f*(X[fluid.IX(1,0)]+X[fluid.IX(0,1)]);
    X[fluid.IX(0,fluid.N+1)]            = 0.5f*(X[fluid.IX(1,fluid.N+1)]+X[fluid.IX(0,fluid.N)]);
    X[fluid.IX(fluid.N+1,0)]            = 0.5f*(X[fluid.IX(fluid.N,0)]+X[fluid.IX(fluid.N+1,1)]);
    X[fluid.IX(fluid.N+1,fluid.N+1)]    = 0.5f*(X[fluid.IX(fluid.N,fluid.N+1)]+X[fluid.IX(fluid.N+1,fluid.N)]);
}

