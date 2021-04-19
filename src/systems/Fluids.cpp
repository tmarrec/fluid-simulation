#include "Fluids.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unistd.h>
#include <omp.h>


void Fluids::init(std::shared_ptr<Renderer> renderer)
{
	_renderer = renderer;
}

void Fluids::update()
{
	for (auto const& entity : mEntities)
	{
		auto& fluid = gCoordinator.GetComponent<Fluid3D>(entity);

		int N = fluid.N/2;

        float p = 256;

	    fluid.substanceField[fluid.IX(N, N, N*2-3)] = p;

        float z = 300;
	    fluid.velocityFieldZ[fluid.IX(N, N, N*2-2)] = -z;
		
		Vstep(fluid);
		Sstep(fluid);
		updateRender(fluid);
	}
}

void Fluids::Vstep(Fluid3D& fluid)
{
	addSource(fluid, fluid.velocityFieldX, fluid.velocityFieldPrevX);
	addSource(fluid, fluid.velocityFieldY, fluid.velocityFieldPrevY);
	addSource(fluid, fluid.velocityFieldZ, fluid.velocityFieldPrevZ);

	swap(fluid.velocityFieldPrevX, fluid.velocityFieldX);
	diffuse(fluid, fluid.velocityFieldX, fluid.velocityFieldPrevX, 1);
	swap(fluid.velocityFieldPrevY, fluid.velocityFieldY);
	diffuse(fluid, fluid.velocityFieldY, fluid.velocityFieldPrevY, 2);
	swap(fluid.velocityFieldPrevZ, fluid.velocityFieldZ);
	diffuse(fluid, fluid.velocityFieldZ, fluid.velocityFieldPrevZ, 3);

	project(fluid, fluid.velocityFieldX, fluid.velocityFieldY, fluid.velocityFieldZ, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY);

	swap(fluid.velocityFieldPrevX, fluid.velocityFieldX);
	swap(fluid.velocityFieldPrevY, fluid.velocityFieldY);
	swap(fluid.velocityFieldPrevZ, fluid.velocityFieldZ);

	advect(fluid, fluid.velocityFieldX, fluid.velocityFieldPrevX, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, 1);
	advect(fluid, fluid.velocityFieldY, fluid.velocityFieldPrevY, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, 2);
	advect(fluid, fluid.velocityFieldZ, fluid.velocityFieldPrevZ, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, 3);

	project(fluid, fluid.velocityFieldX, fluid.velocityFieldY, fluid.velocityFieldZ, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY);
}

void Fluids::Sstep(Fluid3D& fluid)
{
	addSource(fluid, fluid.substanceField, fluid.substanceFieldPrev);
	swap(fluid.substanceFieldPrev, fluid.substanceField);
	diffuse(fluid, fluid.substanceField, fluid.substanceFieldPrev, 0);
	swap(fluid.substanceFieldPrev, fluid.substanceField);
	advect(fluid, fluid.substanceField, fluid.substanceFieldPrev, fluid.velocityFieldX, fluid.velocityFieldY, fluid.velocityFieldZ, 0);
}

void Fluids::addSource(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& S) const
{
	for (std::uint32_t i = 0; i < (fluid.N+2)*(fluid.N+2)*(fluid.N+2); ++i)
	{
		X[i] += fluid.dt * S[i];
	}
}

void Fluids::diffuse(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& Xprev, std::uint8_t b) const
{
	float a = fluid.dt * fluid.viscosity * fluid.N * fluid.N * fluid.N;
	linSolve(fluid, X, Xprev, a, 1+6*a, b);
}

void Fluids::advect(Fluid3D& fluid, std::vector<float>& D, std::vector<float>& Dprev, std::vector<float>& X, std::vector<float>& Y, std::vector<float>& Z, std::uint8_t b) const
{
	float dt0 = fluid.dt * fluid.N;
	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t j = 1; j <= fluid.N; ++j)
		{
			for (std::uint32_t i = 1; i <= fluid.N; ++i)
			{
				float x = i-dt0*X[fluid.IX(i,j,k)];
				float y = j-dt0*Y[fluid.IX(i,j,k)];
				float z = k-dt0*Z[fluid.IX(i,j,k)];
				if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
				{
					ERROR("The simulation time-step is probably too high..");
				}

				x = std::clamp(x, 0.5f, fluid.N + 0.5f);
				std::uint32_t i0 = static_cast<std::uint32_t>(x);
				std::uint32_t i1 = i0 + 1;
				y = std::clamp(y, 0.5f, fluid.N + 0.5f);
				std::uint32_t j0 = static_cast<std::uint32_t>(y);
				std::uint32_t j1 = j0 + 1;
				z = std::clamp(z, 0.5f, fluid.N + 0.5f);
				std::uint32_t k0 = static_cast<std::uint32_t>(z);
				std::uint32_t k1 = k0 + 1;

				float s1 = x	- i0;
				float s0 = 1.0f - s1;
				float t1 = y	- j0;
				float t0 = 1.0f - t1;
				float u1 = z	- k0;
				float u0 = 1.0f - u1;

				D[fluid.IX(i,j,k)] = 
					s0*(t0*(u0*Dprev[fluid.IX(i0,j0,k0)]
							+u1*Dprev[fluid.IX(i0,j0,k1)])
						+(t1*(u0*Dprev[fluid.IX(i0,j1,k0)]
							+u1*Dprev[fluid.IX(i0,j1,k1)])))
					+s1*(t0*(u0*Dprev[fluid.IX(i1,j0,k0)]
							+u1*Dprev[fluid.IX(i1,j0,k1)])
						+(t1*(u0*Dprev[fluid.IX(i1,j1,k0)]
							+u1*Dprev[fluid.IX(i1,j1,k1)])));
			}
		}
	}
	setBnd(fluid, D, b);
}

void Fluids::linSolve(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& Xprev, float a, float c, std::uint8_t b) const
{
	float h = 1.0f/c;
	for (std::uint32_t l = 0; l < fluid.solveNb; l++)
	{
		for (std::uint32_t k = 1; k < fluid.N; k++)
        {
			for (std::uint32_t j = 1; j < fluid.N; j++)
            {
				for (std::uint32_t i = 1; i < fluid.N; i++)
                {
					X[fluid.IX(i,j,k)] =
                        (Xprev[fluid.IX(i,j,k)]
							+a*(X[fluid.IX(i+1,j,k)]+
								X[fluid.IX(i-1,j,k)]+
								X[fluid.IX(i,j+1,k)]+
								X[fluid.IX(i,j-1,k)]+
								X[fluid.IX(i,j,k+1)]+
								X[fluid.IX(i,j,k-1)]
						   ))*h;
				}
			}
		}
		setBnd(fluid, X, b);
	}
}


void Fluids::project(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& Y, std::vector<float>& Z, std::vector<float>& p, std::vector<float>& div) const
{
    float h = 1.0f/fluid.N;
	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t j = 1; j <= fluid.N; ++j)
		{
			for (std::uint32_t i = 1; i <= fluid.N; ++i)
			{
                div[fluid.IX(i,j,k)] = -(1.0f/3.0f)*
                    ((X[fluid.IX(i+1,j,k)]-X[fluid.IX(i-1,j,k)])*h+
                     (Y[fluid.IX(i,j+1,k)]-Y[fluid.IX(i,j-1,k)])*h+
                     (Z[fluid.IX(i,j,k+0)]-Z[fluid.IX(i,j,k-1)])*h);
				p[fluid.IX(i,j,k)] = 0;
			}
		}
	}
	setBnd(fluid, div, 0);
	setBnd(fluid, p, 0);

	linSolve(fluid, p, div, 1, 6, 0);

	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t j = 1; j <= fluid.N; ++j)
		{
			for (std::uint32_t i = 1; i <= fluid.N; ++i)
			{
				X[fluid.IX(i,j,k)] -= 0.5f*fluid.N*(p[fluid.IX(i+1,j,k)]-p[fluid.IX(i-1,j,k)]);
				Y[fluid.IX(i,j,k)] -= 0.5f*fluid.N*(p[fluid.IX(i,j+1,k)]-p[fluid.IX(i,j-1,k)]);
				Z[fluid.IX(i,j,k)] -= 0.5f*fluid.N*(p[fluid.IX(i,j,k+1)]-p[fluid.IX(i,j,k-1)]);
			}
		}
	}

	setBnd(fluid, X, 1);
	setBnd(fluid, Y, 2);
	setBnd(fluid, Z, 3);
}

void Fluids::setBnd(Fluid3D& fluid, std::vector<float>& X, std::uint8_t b) const
{
	// Faces
	for (std::uint32_t j = 1; j <= fluid.N; ++j)
	{
		for (std::uint32_t i = 1; i <= fluid.N; ++i)
		{
			X[fluid.IX(i,j,0)]		    = b == 3 ? -X[fluid.IX(i,j,1)]		    : X[fluid.IX(i,j,1)];
			X[fluid.IX(i,j,fluid.N+1)]  = b == 3 ? -X[fluid.IX(i,j,fluid.N)]	: X[fluid.IX(i,j,fluid.N)];
		}
	}
	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t i = 1; i <= fluid.N; ++i)
		{
			X[fluid.IX(i,0,k)]		    = b == 2 ? -X[fluid.IX(i,1,k)]		    : X[fluid.IX(i,1,k)];
			X[fluid.IX(i,fluid.N+1,k)]  = b == 2 ? -X[fluid.IX(i,fluid.N,k)]	: X[fluid.IX(i,fluid.N,k)];
		}
	}
	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t j = 1; j <= fluid.N; ++j)
		{
			X[fluid.IX(0,j,k)]		    = b == 1 ? -X[fluid.IX(1,j,k)]		    : X[fluid.IX(1,j,k)];
			X[fluid.IX(fluid.N+1,j,k)]  = b == 1 ? -X[fluid.IX(fluid.N,j,k)]	: X[fluid.IX(fluid.N,j,k)];
		}
	}

	// Edges
	for (std::uint32_t i = 1; i <= fluid.N; ++i)
	{
		X[fluid.IX(i,0,0)] = 0.5f*(X[fluid.IX(i,1,0)]+X[fluid.IX(i,0,1)]);
		X[fluid.IX(i,fluid.N+1,0)] = 0.5f*(X[fluid.IX(i,fluid.N,0)]+X[fluid.IX(i,fluid.N+1,1)]);
		X[fluid.IX(i,0,fluid.N+1)] = 0.5f*(X[fluid.IX(i,0,fluid.N)]+X[fluid.IX(i,1,fluid.N+1)]);
		X[fluid.IX(i,fluid.N+1,fluid.N+1)] = 0.5f*(X[fluid.IX(i,fluid.N,fluid.N+1)]+X[fluid.IX(i,fluid.N+1,fluid.N)]);

		X[fluid.IX(0,i,0)] = 0.5f*(X[fluid.IX(1,i,0)]+X[fluid.IX(0,i,1)]);
		X[fluid.IX(fluid.N+1,i,0)] = 0.5f*(X[fluid.IX(fluid.N,i,0)]+X[fluid.IX(fluid.N+1,i,1)]);
		X[fluid.IX(0,i,fluid.N+1)] = 0.5f*(X[fluid.IX(0,i,fluid.N)]+X[fluid.IX(1,i,fluid.N+1)]);
		X[fluid.IX(fluid.N+1,i,fluid.N+1)] = 0.5f*(X[fluid.IX(fluid.N,i,fluid.N+1)]+X[fluid.IX(fluid.N+1,i,0)]);

		X[fluid.IX(0,0,i)] = 0.5f*(X[fluid.IX(0,1,i)]+X[fluid.IX(1,0,i)]);
		X[fluid.IX(0,fluid.N+1,i)] = 0.5f*(X[fluid.IX(0,fluid.N,i)]+X[fluid.IX(1,fluid.N+1,i)]);
		X[fluid.IX(fluid.N+1,0,i)] = 0.5f*(X[fluid.IX(fluid.N,0,i)]+X[fluid.IX(fluid.N+1,1,i)]);
		X[fluid.IX(fluid.N+1,fluid.N+1,i)] = 0.5f*(X[fluid.IX(fluid.N+1,fluid.N,i)]+X[fluid.IX(fluid.N,fluid.N+1,i)]);
	}

	// Corners
	X[fluid.IX(0,0,0)]				    = 0.33f*(X[fluid.IX(1,0,0)]+X[fluid.IX(0,1,0)]+X[fluid.IX(0,0,1)]);
	X[fluid.IX(0,fluid.N+1,0)]		    = 0.33f*(X[fluid.IX(1,fluid.N+1,0)]+X[fluid.IX(0,fluid.N,0)]+X[fluid.IX(0,fluid.N+1,1)]);
	X[fluid.IX(0,0,fluid.N+1)]		    = 0.33f*(X[fluid.IX(1,0,fluid.N+1)]+X[fluid.IX(0,1,fluid.N+1)]+X[fluid.IX(0,0,fluid.N+2)]);
	X[fluid.IX(0,fluid.N+1,fluid.N+1)]  = 0.33f*(X[fluid.IX(1,fluid.N+1,fluid.N+1)]+X[fluid.IX(0,fluid.N,fluid.N+1)]+X[fluid.IX(0,fluid.N+1,fluid.N)]);

	X[fluid.IX(fluid.N+1,0,0)]				    = 0.33f*(X[fluid.IX(fluid.N,0,0)]+X[fluid.IX(fluid.N+1,1,0)]+X[fluid.IX(fluid.N+1,0,1)]);
	X[fluid.IX(fluid.N+1,fluid.N+1,0)]		    = 0.33f*(X[fluid.IX(fluid.N,fluid.N+1,0)]+X[fluid.IX(fluid.N+1,fluid.N,0)]+X[fluid.IX(fluid.N+1,fluid.N+1,1)]);
	X[fluid.IX(fluid.N+1,0,fluid.N+1)]		    = 0.33f*(X[fluid.IX(fluid.N,0,fluid.N+1)]+X[fluid.IX(fluid.N+1,1,fluid.N+1)]+X[fluid.IX(fluid.N+1,0,fluid.N)]);
	X[fluid.IX(fluid.N+1,fluid.N+1,fluid.N+1)]  = 0.33f*(X[fluid.IX(fluid.N,fluid.N+1,fluid.N+1)]+X[fluid.IX(fluid.N+1,fluid.N,fluid.N+1)]+X[fluid.IX(fluid.N+1,fluid.N+1,fluid.N)]);
}

void Fluids::swap(std::vector<float>& X, std::vector<float>& Y) const
{
	X.swap(Y);
}

void Fluids::updateRender(Fluid3D& fluid)
{
	std::vector<std::uint8_t> texture((fluid.N+2)*(fluid.N+2)*(fluid.N+2), 0);

    float maxDensity = 0.0f;
	for (std::uint32_t i = 0; i < (fluid.N+2)*(fluid.N+2)*(fluid.N+2); ++i)
	{
        std::uint32_t density = fluid.substanceField[i];
        if (density > maxDensity)
        {
            maxDensity = density;
        }
	}

	for (std::uint32_t i = 0; i < (fluid.N+2)*(fluid.N+2)*(fluid.N+2); ++i)
	{
		texture[i] = (fluid.substanceField[i]/maxDensity)*255;
	}

	const auto& textureGL = gCoordinator.GetComponent<Material>(fluid.entity).texture;
	_renderer->initTexture3D(texture, textureGL);
}
