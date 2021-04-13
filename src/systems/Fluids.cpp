#include "Fluids.h"

void Fluids::init(std::shared_ptr<Renderer> renderer)
{
    _renderer = renderer;
}

void Fluids::update()
{
    for (auto const& entity : mEntities)
    {
        auto& fluid = gCoordinator.GetComponent<Fluid2D>(entity);
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
        int sum = (int)(fluid.substanceField[i].x + fluid.substanceField[i].y + fluid.substanceField[i].z);
        sum = std::clamp(sum, 0, 255);
        texture[i*3] = sum;
        texture[i*3+1] = sum;
        texture[i*3+2] = sum;
    }
    _renderer->updateTexture(texture);
}

void Fluids::Vstep(Fluid2D& fluid)
{
    addSource(fluid);

    fluid.velocityFieldPrev = fluid.velocityField;

    diffuse(fluid, fluid.velocityField, fluid.velocityFieldPrev);
    project(fluid, fluid.velocityField, fluid.velocityFieldPrev);

    fluid.velocityFieldPrev = fluid.velocityField;

    advect(fluid, fluid.velocityField, fluid.velocityFieldPrev, fluid.velocityFieldPrev);
    project(fluid, fluid.velocityField, fluid.velocityFieldPrev);
}

void Fluids::Sstep(Fluid2D& fluid)
{
    //addSource(fluid);

    fluid.substanceFieldPrev = fluid.substanceField;

    diffuse(fluid, fluid.substanceField, fluid.substanceFieldPrev);

    fluid.substanceFieldPrev = fluid.substanceField;

    advect(fluid, fluid.substanceField, fluid.substanceFieldPrev, fluid.velocityField);
}

void Fluids::addSource(Fluid2D& fluid)
{
    // Testing
    fluid.velocityField[fluid.N/2+fluid.N*1] = glm::vec3{0.0f,1300.0f,0.0f};
    fluid.substanceField[fluid.N/2+fluid.N*2] = fluid.substanceField[fluid.N/2+fluid.N*2] + glm::vec3{0.0f,255.0f,0.0f};
    fluid.substanceField[fluid.N/2+fluid.N*2+1] = fluid.substanceField[fluid.N/2+fluid.N*2+1] + glm::vec3{0.0f,255.0f,0.0f};
}

void Fluids::diffuse(Fluid2D& fluid, std::vector<glm::vec3>& X, std::vector<glm::vec3>& Xprev)
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
        setBnd(fluid, X);
    }
}

void Fluids::advect(Fluid2D& fluid, std::vector<glm::vec3>& D, std::vector<glm::vec3>& Dprev, std::vector<glm::vec3>& U)
{
    float dt0 = fluid.dt * fluid.N;
    for (std::uint32_t i = 1; i <= fluid.N; ++i)
    {
        for (std::uint32_t j = 1; j <= fluid.N; ++j)
        {
            float x = i-dt0*U[fluid.IX(i,j)][0];
            float y = j-dt0*U[fluid.IX(i,j)][1];

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
    setBnd(fluid, D);
}

void Fluids::project(Fluid2D& fluid, std::vector<glm::vec3>& U, std::vector<glm::vec3>& Uprev)
{
    float h = 1.0f/fluid.N;
    for (std::uint32_t i = 1; i <= fluid.N; ++i)
    {
        for (std::uint32_t j = 1; j <= fluid.N; ++j)
        {
            Uprev[fluid.IX(i,j)][1] = -0.5f*h*(U[fluid.IX(i+1,j)][0]-U[fluid.IX(i-1,j)][0]+U[fluid.IX(i,j+1)][1]-U[fluid.IX(i,j-1)][1]);
            Uprev[fluid.IX(i,j)][0] = 0;
        }
    }
    setBnd(fluid, Uprev);

    for (std::uint32_t k = 0; k < 20; ++k)
    {
        for (std::uint32_t i = 1; i <= fluid.N; ++i)
        {
            for (std::uint32_t j = 1; j <= fluid.N; ++j)
            {
                Uprev[fluid.IX(i,j)][0] = (Uprev[fluid.IX(i,j)][1]+Uprev[fluid.IX(i-1,j)][0]+Uprev[fluid.IX(i+1,j)][0]+Uprev[fluid.IX(i,j-1)][0]+Uprev[fluid.IX(i,j+1)][0])/4;
            }
        }
        setBnd(fluid, Uprev);
    }
    
    for (std::uint32_t i = 1; i <= fluid.N; ++i)
    {
        for (std::uint32_t j = 1; j <= fluid.N; ++j)
        {
            U[fluid.IX(i,j)][0] -= 0.5f*(Uprev[fluid.IX(i+1,j)][0]-Uprev[fluid.IX(i-1,j)][0])/h;
            U[fluid.IX(i,j)][1] -= 0.5f*(Uprev[fluid.IX(i,j+1)][0]-Uprev[fluid.IX(i,j-1)][0])/h;
        }
    }
    setBnd(fluid, U);
}

void Fluids::setBnd(Fluid2D& fluid, std::vector<glm::vec3>& U)
{
    for (std::uint32_t i = 1; i <= fluid.N; ++i)
    {
        U[fluid.IX(0,i)]            = glm::vec3{0,0,0};
        U[fluid.IX(i,0)]            = glm::vec3{0,0,0};
        U[fluid.IX(fluid.N+1,i)]    = glm::vec3{0,0,0};
        U[fluid.IX(i,fluid.N+1)]    = glm::vec3{0,0,0};
    }
    U[fluid.IX(0,0)]                    = 0.5f*(U[fluid.IX(1,0)]+U[fluid.IX(0,1)]);
    U[fluid.IX(0,fluid.N+1)]            = 0.5f*(U[fluid.IX(1,fluid.N+1)]+U[fluid.IX(0,fluid.N)]);
    U[fluid.IX(fluid.N+1,0)]            = 0.5f*(U[fluid.IX(fluid.N,0)]+U[fluid.IX(fluid.N+1,1)]);
    U[fluid.IX(fluid.N+1,fluid.N+1)]    = 0.5f*(U[fluid.IX(fluid.N,fluid.N+1)]+U[fluid.IX(fluid.N+1,fluid.N)]);
}
