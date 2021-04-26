#include "Fluids.h"

#include <Eigen/Sparse>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/IterativeLinearSolvers/BasicPreconditioners.h>
#include <cstdint>

void Fluids::init(std::shared_ptr<Renderer> renderer)
{
	_renderer = renderer;
    std::srand(std::time(nullptr));
}


void Fluids::reset(bool force)
{
#ifdef DEBUG_GUI
	for (auto const& entity : mEntities)
	{
		auto& fluid = gCoordinator.GetComponent<Fluid3D>(entity);
        bool rebuild = _debugDt != fluid.dt || _debugViscosity != fluid.viscosity || _debugN != (int)fluid.N || fluid.diffusion != _debugDiffusion;
        if (force || rebuild)
        {
            fluid = Fluid3D();
            fluid.entity = entity;
            fluid.N = _debugN;
            fluid.viscosity = _debugViscosity;
            fluid.diffusion = _debugDiffusion;
            fluid.dt = _debugDt;
            fluid.init();
            _initCG = false;
        }
		auto& material = gCoordinator.GetComponent<Material>(entity);
        material.absorption = _debugAbsorption;
        material.lightIntensity.x = _debugLightIntensity[0];
        material.lightIntensity.y = _debugLightIntensity[1];
        material.lightIntensity.z = _debugLightIntensity[2];
    }
#endif
}

void Fluids::update()
{
    if (mEntities.size() > 1)
    {
        ERROR("Multiple fluids not supported.");
    }
	for (auto const& entity : mEntities)
	{
		auto& fluid = gCoordinator.GetComponent<Fluid3D>(entity);
        
        if (!_initCG)
        {
            _cgProject.compute(fluid.laplacianProject);
            _cgDiffuse.compute(fluid.laplacianDiffuse);
            _cgViscosity.compute(fluid.laplacianViscosity);
            _initCG = true;
        }

        /* Testings */
		int N = fluid.N/2;

        float p = 8;
	    fluid.substanceField[fluid.IX(N, N, N*2-2)] = p;
	    fluid.substanceField[fluid.IX(N, N, 2)] = p;

        float z = 512;
	    fluid.velocityFieldZ[fluid.IX(N, N, N*2-2)] = -z;
	    fluid.velocityFieldZ[fluid.IX(N, N, 2)] = z;
        /* End Testings */

#ifdef DEBUG_GUI
        auto VstepTimeStart = std::chrono::high_resolution_clock::now();
#endif
		Vstep(fluid);
#ifdef DEBUG_GUI
        auto VstepTimeStop = std::chrono::high_resolution_clock::now();
        _debugVstepTimes.emplace_back(std::chrono::duration<float, std::chrono::seconds::period>(VstepTimeStop - VstepTimeStart).count());
        if (_debugVstepTimes.size() > TIME_ECHANT_NB)
            _debugVstepTimes.pop_front();
        auto SstepTimeStart = std::chrono::high_resolution_clock::now();
#endif
		Sstep(fluid);
#ifdef DEBUG_GUI
        auto SstepTimeStop = std::chrono::high_resolution_clock::now();
        _debugSstepTimes.emplace_back(std::chrono::duration<float, std::chrono::seconds::period>(SstepTimeStop - SstepTimeStart).count());
        if (_debugSstepTimes.size() > TIME_ECHANT_NB)
            _debugSstepTimes.pop_front();
        auto TextureTimeStart = std::chrono::high_resolution_clock::now();
#endif
		updateRender(fluid);
#ifdef DEBUG_GUI
        auto TextureTimeStop = std::chrono::high_resolution_clock::now();
        _debugTextureTimes.emplace_back(std::chrono::duration<float, std::chrono::seconds::period>(TextureTimeStop - TextureTimeStart).count());
        if (_debugTextureTimes.size() > TIME_ECHANT_NB)
            _debugTextureTimes.pop_front();
#endif
	}
}

void Fluids::Vstep(Fluid3D& fluid)
{
    /*
	addSource(fluid, fluid.velocityFieldX, fluid.velocityFieldPrevX);
	addSource(fluid, fluid.velocityFieldY, fluid.velocityFieldPrevY);
	addSource(fluid, fluid.velocityFieldZ, fluid.velocityFieldPrevZ);
    */

#ifdef DEBUG_GUI
    auto VstepDiffuseTimeStart = std::chrono::high_resolution_clock::now();
#endif
	diffuse(fluid, fluid.velocityFieldPrevX, fluid.velocityFieldX, 1, _cgViscosity);
	diffuse(fluid, fluid.velocityFieldPrevY, fluid.velocityFieldY, 2, _cgViscosity);
	diffuse(fluid, fluid.velocityFieldPrevZ, fluid.velocityFieldZ, 3, _cgViscosity);
#ifdef DEBUG_GUI
    auto VstepDiffuseTimeStop = std::chrono::high_resolution_clock::now();
    _debugVstepDiffuseTimes.emplace_back(std::chrono::duration<float, std::chrono::seconds::period>(VstepDiffuseTimeStop - VstepDiffuseTimeStart).count());
    if (_debugVstepDiffuseTimes.size() > TIME_ECHANT_NB)
        _debugVstepDiffuseTimes.pop_front();

    auto VstepProjectTimeStart = std::chrono::high_resolution_clock::now();
#endif
	project(fluid, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, fluid.velocityFieldX, fluid.velocityFieldY);
#ifdef DEBUG_GUI
    auto VstepProjectTimeStop = std::chrono::high_resolution_clock::now();
    _debugVstepProjectTimes.emplace_back(std::chrono::duration<float, std::chrono::seconds::period>(VstepProjectTimeStop - VstepProjectTimeStart).count());
    if (_debugVstepProjectTimes.size() > TIME_ECHANT_NB)
        _debugVstepProjectTimes.pop_front();

    auto VstepAdvectTimeStart = std::chrono::high_resolution_clock::now();
#endif
	advect(fluid, fluid.velocityFieldX, fluid.velocityFieldPrevX, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, 1);
	advect(fluid, fluid.velocityFieldY, fluid.velocityFieldPrevY, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, 2);
	advect(fluid, fluid.velocityFieldZ, fluid.velocityFieldPrevZ, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, 3);
#ifdef DEBUG_GUI
    auto VstepAdvectTimeStop = std::chrono::high_resolution_clock::now();
    _debugVstepAdvectTimes.emplace_back(std::chrono::duration<float, std::chrono::seconds::period>(VstepAdvectTimeStop - VstepAdvectTimeStart).count());
    if (_debugVstepAdvectTimes.size() > TIME_ECHANT_NB)
        _debugVstepAdvectTimes.pop_front();

    VstepProjectTimeStart = std::chrono::high_resolution_clock::now();
#endif
	project(fluid, fluid.velocityFieldX, fluid.velocityFieldY, fluid.velocityFieldZ, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY);
#ifdef DEBUG_GUI
    VstepProjectTimeStop = std::chrono::high_resolution_clock::now();
    _debugVstepProjectTimes.back() += std::chrono::duration<float, std::chrono::seconds::period>(VstepProjectTimeStop - VstepProjectTimeStart).count();
#endif
}

void Fluids::Sstep(Fluid3D& fluid)
{
	//addSource(fluid, fluid.substanceField, fluid.substanceFieldPrev);
    
#ifdef DEBUG_GUI
    auto SstepDiffuseTimeStart = std::chrono::high_resolution_clock::now();
#endif
	diffuse(fluid, fluid.substanceFieldPrev, fluid.substanceField, 0, _cgDiffuse);

#ifdef DEBUG_GUI
    auto SstepDiffuseTimeStop = std::chrono::high_resolution_clock::now();
    _debugSstepDiffuseTimes.emplace_back(std::chrono::duration<float, std::chrono::seconds::period>(SstepDiffuseTimeStop - SstepDiffuseTimeStart).count());
    if (_debugSstepDiffuseTimes.size() > TIME_ECHANT_NB)
        _debugSstepDiffuseTimes.pop_front();

    auto SstepAdvectTimeStart = std::chrono::high_resolution_clock::now();
#endif
	advect(fluid, fluid.substanceField, fluid.substanceFieldPrev, fluid.velocityFieldX, fluid.velocityFieldY, fluid.velocityFieldZ, 0);
#ifdef DEBUG_GUI
    auto SstepAdvectTimeStop = std::chrono::high_resolution_clock::now();
    _debugSstepAdvectTimes.emplace_back(std::chrono::duration<float, std::chrono::seconds::period>(SstepAdvectTimeStop - SstepAdvectTimeStart).count());
    if (_debugSstepAdvectTimes.size() > TIME_ECHANT_NB)
        _debugSstepAdvectTimes.pop_front();
#endif
}

void Fluids::addSource(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& S) const
{
	for (std::uint32_t i = 0; i < (fluid.N+2)*(fluid.N+2)*(fluid.N+2); ++i)
	{
		X[i] += fluid.dt * S[i];
	}
}

void Fluids::diffuse(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& Xprev, std::uint8_t b, Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper>& cg)
{
    ConjugateGradientMethodLinSolve(fluid, X, Xprev, b, cg);
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

void Fluids::ConjugateGradientMethodLinSolve(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& Xprev, std::uint8_t bs, Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper>& cg)
{
    std::uint32_t n = (fluid.N)*(fluid.N)*(fluid.N);
    Eigen::VectorXd x(n);
    Eigen::VectorXd b(n);

    // Filling matrices
    std::uint32_t bIt = 0;
    for (std::uint32_t k = 1; k <= fluid.N; k++)
    {
        for (std::uint32_t j = 1; j <= fluid.N; j++)
        {
            for (std::uint32_t i = 1; i <= fluid.N; i++)
            {
                b.coeffRef(bIt) = -(Xprev[fluid.IX(i,j,k)]);
                bIt++;
            }
        }
    }
    
    // Solve the linear system
    x = cg.solve(b);
    
    // Write the results
    for (std::uint32_t k = 1; k <= fluid.N; k++)
    {
        for (std::uint32_t j = 1; j <= fluid.N; j++)
        {
            for (std::uint32_t i = 1; i <= fluid.N; i++)
            {
                X[fluid.IX(i,j,k)] = x.coeffRef((i-1)+(j-1)*fluid.N+(k-1)*fluid.N*fluid.N);
            }
        }
    }

    setBnd(fluid, X, bs);
}

void Fluids::project(Fluid3D& fluid, std::vector<float>& X, std::vector<float>& Y, std::vector<float>& Z, std::vector<float>& p, std::vector<float>& div)
{
    float h = 1.0f/fluid.N;
	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t j = 1; j <= fluid.N; ++j)
		{
			for (std::uint32_t i = 1; i <= fluid.N; ++i)
			{
                div[fluid.IX(i,j,k)] = -(0.33f)*
                    ((X[fluid.IX(i+1,j,k)]-X[fluid.IX(i-1,j,k)])*h+
                     (Y[fluid.IX(i,j+1,k)]-Y[fluid.IX(i,j-1,k)])*h+
                     (Z[fluid.IX(i,j,k+1)]-Z[fluid.IX(i,j,k-1)])*h);
				p[fluid.IX(i,j,k)] = 0;
			}
		}
	}
	setBnd(fluid, div, 0);
	setBnd(fluid, p, 0);

    ConjugateGradientMethodLinSolve(fluid, p, div, 0, _cgProject);

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
			X[fluid.IX(i,j,fluid.N)]  = b == 3 ? -X[fluid.IX(i,j,fluid.N-1)]	: X[fluid.IX(i,j,fluid.N-1)];
		}
	}
	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t i = 1; i <= fluid.N; ++i)
		{
			X[fluid.IX(i,0,k)]		    = b == 2 ? -X[fluid.IX(i,1,k)]		    : X[fluid.IX(i,1,k)];
			X[fluid.IX(i,fluid.N,k)]  = b == 2 ? -X[fluid.IX(i,fluid.N-1,k)]	: X[fluid.IX(i,fluid.N-1,k)];
		}
	}
	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t j = 1; j <= fluid.N; ++j)
		{
			X[fluid.IX(0,j,k)]		    = b == 1 ? -X[fluid.IX(1,j,k)]		    : X[fluid.IX(1,j,k)];
			X[fluid.IX(fluid.N,j,k)]  = b == 1 ? -X[fluid.IX(fluid.N-1,j,k)]	: X[fluid.IX(fluid.N-1,j,k)];
		}
	}

	// Edges
	for (std::uint32_t i = 1; i <= fluid.N; ++i)
	{
		X[fluid.IX(i,0,0)]              = 0.5f*(X[fluid.IX(i,1,0)]              +X[fluid.IX(i,0,1)]);
		X[fluid.IX(i,fluid.N,0)]        = 0.5f*(X[fluid.IX(i,fluid.N-1,0)]      +X[fluid.IX(i,fluid.N,1)]);
		X[fluid.IX(i,0,fluid.N)]        = 0.5f*(X[fluid.IX(i,0,fluid.N-1)]      +X[fluid.IX(i,1,fluid.N)]);
		X[fluid.IX(i,fluid.N,fluid.N)]  = 0.5f*(X[fluid.IX(i,fluid.N-1,fluid.N)]+X[fluid.IX(i,fluid.N,fluid.N-1)]);

		X[fluid.IX(0,i,0)]              = 0.5f*(X[fluid.IX(1,i,0)]              +X[fluid.IX(0,i,1)]);
		X[fluid.IX(fluid.N,i,0)]        = 0.5f*(X[fluid.IX(fluid.N-1,i,0)]      +X[fluid.IX(fluid.N,i,1)]);
		X[fluid.IX(0,i,fluid.N)]        = 0.5f*(X[fluid.IX(0,i,fluid.N-1)]      +X[fluid.IX(1,i,fluid.N)]);
		X[fluid.IX(fluid.N,i,fluid.N)]  = 0.5f*(X[fluid.IX(fluid.N-1,i,fluid.N)]+X[fluid.IX(fluid.N,i,0)]);

		X[fluid.IX(0,0,i)]              = 0.5f*(X[fluid.IX(0,1,i)]              +X[fluid.IX(1,0,i)]);
		X[fluid.IX(0,fluid.N,i)]        = 0.5f*(X[fluid.IX(0,fluid.N-1,i)]      +X[fluid.IX(1,fluid.N,i)]);
		X[fluid.IX(fluid.N,0,i)]        = 0.5f*(X[fluid.IX(fluid.N-1,0,i)]      +X[fluid.IX(fluid.N,1,i)]);
		X[fluid.IX(fluid.N,fluid.N,i)]  = 0.5f*(X[fluid.IX(fluid.N,fluid.N-1,i)]+X[fluid.IX(fluid.N-1,fluid.N,i)]);
	}

	// Corners
	X[fluid.IX(0,0,0)]				    = 0.33f*(X[fluid.IX(1,0,0)]             +X[fluid.IX(0,1,0)]                 +X[fluid.IX(0,0,1)]);
	X[fluid.IX(0,fluid.N,0)]		    = 0.33f*(X[fluid.IX(1,fluid.N,0)]       +X[fluid.IX(0,fluid.N-1,0)]         +X[fluid.IX(0,fluid.N,1)]);
	X[fluid.IX(0,0,fluid.N)]		    = 0.33f*(X[fluid.IX(1,0,fluid.N)]       +X[fluid.IX(0,1,fluid.N)]           +X[fluid.IX(0,0,fluid.N+1)]);
	X[fluid.IX(0,fluid.N,fluid.N)]      = 0.33f*(X[fluid.IX(1,fluid.N,fluid.N)] +X[fluid.IX(0,fluid.N-1,fluid.N)]   +X[fluid.IX(0,fluid.N,fluid.N-1)]);

	X[fluid.IX(fluid.N,0,0)]				= 0.33f*(X[fluid.IX(fluid.N-1,0,0)]             +X[fluid.IX(fluid.N,1,0)]               +X[fluid.IX(fluid.N,0,1)]);
	X[fluid.IX(fluid.N,fluid.N,0)]		    = 0.33f*(X[fluid.IX(fluid.N-1,fluid.N,0)]       +X[fluid.IX(fluid.N,fluid.N-1,0)]       +X[fluid.IX(fluid.N,fluid.N,1)]);
	X[fluid.IX(fluid.N,0,fluid.N)]		    = 0.33f*(X[fluid.IX(fluid.N-1,0,fluid.N)]       +X[fluid.IX(fluid.N,1,fluid.N)]         +X[fluid.IX(fluid.N,0,fluid.N-1)]);
	X[fluid.IX(fluid.N,fluid.N,fluid.N)]    = 0.33f*(X[fluid.IX(fluid.N-1,fluid.N,fluid.N)] +X[fluid.IX(fluid.N,fluid.N-1,fluid.N)] +X[fluid.IX(fluid.N,fluid.N,fluid.N-1)]);
}

void Fluids::updateRender(Fluid3D& fluid)
{
	std::vector<std::uint8_t> texture((fluid.N+2)*(fluid.N+2)*(fluid.N+2), 0);

    float maxDensity = 0.0f;
    float sumDensity = 0.0f;
	for (std::uint32_t i = 0; i < (fluid.N+2)*(fluid.N+2)*(fluid.N); ++i)
	{
        std::uint32_t density = fluid.substanceField[i];
        if (density > maxDensity)
        {
            maxDensity = density;
        }
        sumDensity += density;
		//texture[i] = std::clamp(fluid.substanceField[i], 0.0f, 255.0f);;
	}
    //std::cout << sumDensity << std::endl;

	for (std::uint32_t i = 0; i < (fluid.N+2)*(fluid.N+2)*(fluid.N+2); ++i)
	{
		texture[i] = (fluid.substanceField[i]/maxDensity)*255;
	}

	const auto& textureGL = gCoordinator.GetComponent<Material>(fluid.entity).texture;
	_renderer->initTexture3D(texture, textureGL);
}


#ifdef DEBUG_GUI
void Fluids::fluidSetupDebug()
{
	for (auto const& entity : mEntities)
    {
		auto& fluid = gCoordinator.GetComponent<Fluid3D>(entity);
        _debugViscosity = fluid.viscosity;
        _debugDiffusion = fluid.diffusion;
        _debugDt = fluid.dt;
        _debugN = fluid.N;

		auto& material = gCoordinator.GetComponent<Material>(entity);
        _debugLightIntensity[0] = material.lightIntensity.x;
        _debugLightIntensity[1] = material.lightIntensity.y;
        _debugLightIntensity[2] = material.lightIntensity.z;
        _debugAbsorption = material.absorption;
    }
}

void Fluids::fluidDebugTool()
{
    ImGui::Begin("Fluid");
    ImGui::SliderInt("N", &_debugN, 4, 40);
    ImGui::SliderFloat("viscosity", &_debugViscosity, 0.0f, 16.0f);
    ImGui::SliderFloat("diffusion", &_debugDiffusion, 0.0f, 16.0f);
    ImGui::SliderFloat("dt", &_debugDt, 0.000001f, 0.0002f, "%.8f");
    ImGui::SliderFloat("absorption", &_debugAbsorption, 1.0f, 100.0f);
    ImGui::SliderFloat3("lightIntensity", _debugLightIntensity, 0.0f, 1.0f);
    if (ImGui::Button("Apply"))
    {
        reset();
    }
    if (ImGui::Button("Reset"))
    {
        reset(true);
    }

    ImPlot::SetNextPlotLimits(-0.5f, 0.5f, -0.5f, 0.5f);
    if (ImPlot::BeginPlot("Simulation time", NULL, NULL, ImVec2(175,175)))
    {
        const char* labels[] = {"sstep","vstep","texture"};
        float meanSstepTime = 0.0f;
        for (const auto& t : _debugSstepTimes)
            meanSstepTime += t;
        meanSstepTime /= TIME_ECHANT_NB;
        float meanVstepTime = 0.0f;
        for (const auto& t : _debugVstepTimes)
            meanVstepTime += t;
        meanVstepTime /= TIME_ECHANT_NB;
        float meanTextureTime = 0.0f;
        for (const auto& t : _debugTextureTimes)
            meanTextureTime += t;
        meanTextureTime /= TIME_ECHANT_NB;
        float data[] = {meanSstepTime,meanVstepTime,meanTextureTime};
        ImPlot::PlotPieChart(labels, data, 3, 0.0f, 0.0f, 0.5f, true, "");
        ImPlot::EndPlot();
    }

    ImPlot::SetNextPlotLimits(-0.5f, 0.5f, -0.5f, 0.5f);
    if (ImPlot::BeginPlot("vstep time", NULL, NULL, ImVec2(175,175)))
    {
        const char* labels[] = {"diffuse","advect","project"};
        float meanVstepDiffuseTime = 0.0f;
        for (const auto& t : _debugVstepDiffuseTimes)
            meanVstepDiffuseTime += t;
        meanVstepDiffuseTime /= TIME_ECHANT_NB;
        float meanVstepAdvectTime = 0.0f;
        for (const auto& t : _debugVstepAdvectTimes)
            meanVstepAdvectTime += t;
        meanVstepAdvectTime /= TIME_ECHANT_NB;
        float meanVstepProjectTime = 0.0f;
        for (const auto& t : _debugVstepProjectTimes)
            meanVstepProjectTime += t;
        meanVstepProjectTime /= TIME_ECHANT_NB;
        float data[] = {meanVstepDiffuseTime,meanVstepAdvectTime,meanVstepProjectTime};
        ImPlot::PlotPieChart(labels, data, 3, 0.0f, 0.0f, 0.5f, true, "");
        ImPlot::EndPlot();
    }

    ImPlot::SetNextPlotLimits(-0.5f, 0.5f, -0.5f, 0.5f);
    if (ImPlot::BeginPlot("sstep time", NULL, NULL, ImVec2(175,175)))
    {
        const char* labels[] = {"diffuse","advect"};
        float meanSstepDiffuseTime = 0.0f;
        for (const auto& t : _debugSstepDiffuseTimes)
            meanSstepDiffuseTime += t;
        meanSstepDiffuseTime /= TIME_ECHANT_NB;
        float meanSstepAdvectTime = 0.0f;
        for (const auto& t : _debugSstepAdvectTimes)
            meanSstepAdvectTime += t;
        meanSstepAdvectTime /= TIME_ECHANT_NB;
        float data[] = {meanSstepDiffuseTime,meanSstepAdvectTime};
        ImPlot::PlotPieChart(labels, data, 2, 0.0f, 0.0f, 0.5f, true, "");
        ImPlot::EndPlot();
    }

    ImGui::End();
}
#endif
