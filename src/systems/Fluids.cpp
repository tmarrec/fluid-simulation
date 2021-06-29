#include "Fluids.h"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <Eigen/src/SparseLU/SparseLU.h>
#include <cmath>
#include <cstdint>
#include <iomanip>

void Fluids::init(std::shared_ptr<Renderer> renderer)
{
	_renderer = renderer;
    std::srand(std::time(nullptr));
}

void Fluids::update([[maybe_unused]] std::uint64_t iteration)
{
    if (mEntities.size() > 1)
    {
        ERROR("Multiple fluids not supported.");
    }
	for (auto const& entity : mEntities)
	{
		auto& fluid = gCoordinator.GetComponent<Fluid3D>(entity);
        
        /* Testings */
		int N = fluid.N/2;

        float p = 256;
        float z = 512;

        for (std::uint16_t j = 0; j <= fluid.N+2; ++j)
        {
            for (std::uint16_t i = 0; i <= fluid.N+2; ++i)
            {
                //fluid.velocityFieldX[fluid.IX(i, j, 0, 1)] = -50;
                //fluid.velocityFieldY[fluid.IX(i, j, 0, 2)] = 50;
                /*
                std::uint64_t shift = (fluid.N+2)/2;
                double x = (double)i;
                double y = (double)j;
                fluid.velocityFieldX[fluid.IX(i, j, 0)] = (y-shift)*5;
                fluid.velocityFieldY[fluid.IX(i, j, 0)] = -(x-shift)*5;
                */
            }
        }

        fluid.substanceField[fluid.IX(N-(N/2), N-(N/2), 0, 0)] = p;
        fluid.velocityFieldX[fluid.IX(N-(N/2), N-(N/2), 0, 1)] = z;
        fluid.velocityFieldY[fluid.IX(N-(N/2), N-(N/2), 0, 2)] = z;

        fluid.substanceField[fluid.IX(N+(N/2), N+(N/2), 0, 0)] = p;
        fluid.velocityFieldX[fluid.IX(N+(N/2), N+(N/2), 0, 1)] = -z;
        fluid.velocityFieldY[fluid.IX(N+(N/2), N+(N/2), 0, 2)] = -z;

        fluid.substanceField[fluid.IX(N+(N/2), N-(N/2), 0, 0)] = p;
        fluid.velocityFieldX[fluid.IX(N+(N/2), N-(N/2), 0, 1)] = -z;
        fluid.velocityFieldY[fluid.IX(N+(N/2), N-(N/2), 0, 2)] = z;

        fluid.substanceField[fluid.IX(N-(N/2), N+(N/2), 0, 0)] = p;
        fluid.velocityFieldX[fluid.IX(N-(N/2), N+(N/2), 0, 1)] = z;
        fluid.velocityFieldY[fluid.IX(N-(N/2), N+(N/2), 0, 2)] = -z;


        /* End Testings */

        auto start = std::chrono::high_resolution_clock::now();
        Vstep(fluid);
        auto end = std::chrono::high_resolution_clock::now();
        VstepTime += std::chrono::duration<float, std::chrono::seconds::period>(end - start).count();
        start = std::chrono::high_resolution_clock::now();
        Sstep(fluid);
        end = std::chrono::high_resolution_clock::now();
        SstepTime += std::chrono::duration<float, std::chrono::seconds::period>(end - start).count();
        updateRender(fluid);

        _renderer->updateDynamicLine(fluid.velocityFieldX, fluid.velocityFieldY);

        /*
        std::cout << std::fixed << std::setprecision(2) << std::endl;
        std::cout << std::endl << "=== X ===" << std::endl;
        for (std::uint64_t j = 0; j < fluid.N+2; ++j)
        {
            for (std::uint64_t i = 0; i < fluid.N+2+1; ++i)
            {
                if (fluid.velocityFieldX[fluid.IX(i,j,0,1)] >= 0)
                    std::cout << " ";
                std::cout << fluid.velocityFieldX[fluid.IX(i,j,0,1)] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << "=== Y ===" << std::endl;
        for (std::uint64_t j = 0; j < fluid.N+2+1; ++j)
        {
            for (std::uint64_t i = 0; i < fluid.N+2; ++i)
            {
                if (fluid.velocityFieldY[fluid.IX(i,j,0,2)] >= 0)
                    std::cout << " ";
                std::cout << fluid.velocityFieldY[fluid.IX(i,j,0,2)] << " ";
            }
            std::cout << std::endl;
        }
        */


        // LEVEL-SET
        /*
        double r = 10;
        if (iteration == 0)
        {
            glm::vec2 p;
            p.x = 42;
            p.y = 42;
            particles.emplace_back(p);
        }
        
        if (iteration == 0)
        {
            for (std::uint64_t j = 0; j < fluid.N+2; ++j)
            {
                for (std::uint64_t i = 0; i < fluid.N+2; ++i)
                {
                    double dist = std::sqrt(std::pow(i-particles.front().x,2)+std::pow(j-particles.front().y,2))-r;
                    for (const auto& p : particles)
                    {
                        dist = std::min(dist, std::sqrt(std::pow(i-p.x,2)+std::pow(j-p.y,2))-r);
                        dist = cos((float)i/4)+cos((float)j/4);
                    }
                    fluid.implicitFunctionFieldPrev[fluid.IX(i,j,0)] = dist;
                }
            }
        }
        advect(fluid, fluid.implicitFunctionField, fluid.implicitFunctionFieldPrev, fluid.velocityFieldX, fluid.velocityFieldY, fluid.velocityFieldZ, 0);
        fluid.implicitFunctionFieldPrev = fluid.implicitFunctionField;

        reinitLevelSet(fluid, 64);
        */


        /*
        std::cout << "Vstep: " << VstepTime/(iteration+1) << ", Sstep: " << SstepTime/(iteration+1) << std::endl;
        std::cout << "VstepDiffuse: " << VstepDiffuseTime/(iteration+1) << ", VstepProject: " << VstepProjectTime/(iteration+1) << ", VstepAdvect: " << VstepAdvectTime/(iteration+1) << std::endl;
        std::cout << "SstepDiffuse: " << SstepDiffuseTime/(iteration+1) << ", SstepAdvect: " << SstepAdvectTime/(iteration+1) << std::endl;
        */
	}
}

void Fluids::reinitLevelSet(Fluid3D& fluid, const std::uint64_t nbIte) const
{
    /*
    std::vector<double> Ssf = fluid.implicitFunctionFieldPrev;
    std::vector<double> n = fluid.implicitFunctionFieldPrev;
    double dx = 1.0/fluid.N;

    // Init smoothing function
    for (std::uint64_t j = 0; j < fluid.N+2; ++j)
    {
        for (std::uint64_t i = 0; i <= fluid.N+2; ++i)
        {
            double O0 = fluid.implicitFunctionFieldPrev[fluid.IX(i,j,0)];
            Ssf[fluid.IX(i,j,0)] = O0 / (std::sqrt(std::pow(O0, 2) + std::pow(1.0/fluid.N, 2)));
        }
    }

    // Step forward in fictious time
    for (std::uint64_t relaxit = 0; relaxit < nbIte; ++relaxit)
    {
        for (std::uint64_t j = 0; j < fluid.N+2; ++j)
        {
            for (std::uint64_t i = 0; i < fluid.N+2; ++i)
            {
                double gO = gradLength(fluid, fluid.implicitFunctionFieldPrev, i, j);
                n[fluid.IX(i,j,0)] = fluid.implicitFunctionFieldPrev[fluid.IX(i,j,0)] + (0.5 * dx * (- Ssf[fluid.IX(i,j,0)] * (gO - 1.0)));
            }
        }
        fluid.implicitFunctionFieldPrev = n;
    }
    */
}

double Fluids::gradLength(const Fluid3D& fluid, const std::vector<double>& X, const std::uint64_t i, const std::uint64_t j) const
{
    /*
    double gradI = 0;
    double gradJ = 0;
    if (i == 0)
    {
        gradI = X[fluid.IX(1,j,0)] - X[fluid.IX(0,j,0)];
    }
    else if (i == fluid.N+1)
    {
        gradI = X[fluid.IX(fluid.N+1,j,0)] - X[fluid.IX(fluid.N,j,0)];
    }
    else
    {
        if (std::abs(X[fluid.IX(i+1,j,0)]) < std::abs(X[fluid.IX(i-1,j,0)]))
        {
            gradI = X[fluid.IX(i,j,0)] - X[fluid.IX(i+1,j,0)];
        }
        else
        {
            gradI = X[fluid.IX(i-1,j,0)] - X[fluid.IX(i,j,0)];
        }
    }
    if (j == 0)
    {
        gradJ = X[fluid.IX(i,1,0)] - X[fluid.IX(i,0,0)];
    }
    else if (j == fluid.N+1)
    {
        gradJ = X[fluid.IX(i,fluid.N+1,0)] - X[fluid.IX(i,fluid.N,0)];
    }
    else
    {
        if (std::abs(X[fluid.IX(i,j+1,0)]) < std::abs(X[fluid.IX(i,j-1,0)]))
        {
            gradJ = X[fluid.IX(i,j,0)] - X[fluid.IX(i,j+1,0)];
        }
        else
        {
            gradJ = X[fluid.IX(i,j-1,0)] - X[fluid.IX(i,j,0)];
        }
    }
    return std::sqrt(std::pow(gradI, 2)+std::pow(gradJ, 2));
    */
}

void Fluids::Vstep(Fluid3D& fluid)
{
    /*
	addSource(fluid, fluid.velocityFieldX, fluid.velocityFieldPrevX);
	addSource(fluid, fluid.velocityFieldY, fluid.velocityFieldPrevY);
	addSource(fluid, fluid.velocityFieldZ, fluid.velocityFieldPrevZ);
    */

    auto start = std::chrono::high_resolution_clock::now();
    diffuse(fluid, fluid.velocityFieldPrevX, fluid.velocityFieldX, 1, fluid.laplacianViscosity);
    diffuse(fluid, fluid.velocityFieldPrevY, fluid.velocityFieldY, 2, fluid.laplacianViscosity);
    auto end = std::chrono::high_resolution_clock::now();
    VstepDiffuseTime += std::chrono::duration<float, std::chrono::seconds::period>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    project(fluid, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, fluid.velocityFieldX, fluid.velocityFieldY);
    end = std::chrono::high_resolution_clock::now();
    VstepProjectTime += std::chrono::duration<float, std::chrono::seconds::period>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    advect(fluid, fluid.velocityFieldX, fluid.velocityFieldPrevX, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, 1);
    advect(fluid, fluid.velocityFieldY, fluid.velocityFieldPrevY, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, 2);
    end = std::chrono::high_resolution_clock::now();
    VstepAdvectTime += std::chrono::duration<float, std::chrono::seconds::period>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    project(fluid, fluid.velocityFieldX, fluid.velocityFieldY, fluid.velocityFieldZ, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY);
    end = std::chrono::high_resolution_clock::now();
    VstepProjectTime += std::chrono::duration<float, std::chrono::seconds::period>(end - start).count();
}

void Fluids::Sstep(Fluid3D& fluid)
{
	//addSource(fluid, fluid.substanceField, fluid.substanceFieldPrev);
    auto start = std::chrono::high_resolution_clock::now();
    diffuse(fluid, fluid.substanceFieldPrev, fluid.substanceField, 0, fluid.laplacianDiffuse);
    auto end = std::chrono::high_resolution_clock::now();
    SstepDiffuseTime += std::chrono::duration<float, std::chrono::seconds::period>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    advect(fluid, fluid.substanceField, fluid.substanceFieldPrev, fluid.velocityFieldX, fluid.velocityFieldY, fluid.velocityFieldZ, 0);
    end = std::chrono::high_resolution_clock::now();
    SstepAdvectTime += std::chrono::duration<float, std::chrono::seconds::period>(end - start).count();
}

void Fluids::addSource(const Fluid3D& fluid, std::vector<double>& X, const std::vector<double>& S) const
{
    const std::uint64_t N = fluid.N; 
    const std::uint64_t N32 = (N+2)*(N+2)*(N+2); 
	for (std::uint64_t i = 0; i < N32; ++i)
	{
		X[i] += fluid.dt * S[i];
	}
}

void Fluids::diffuse(const Fluid3D& fluid, std::vector<double>& X, const std::vector<double>& Xprev, const std::uint8_t b, const Laplacian& A)
{
    if (fluid.solver == GAUSS_SEIDEL)
    {
        const double a = b == 0 ? fluid.dt * fluid.diffusion * fluid.N : fluid.dt * fluid.viscosity * fluid.N;
        GaussSeidelRelaxationLinSolve(fluid, X, Xprev, a, fluid.is2D ? 1+4*a : 1+6*a, b);
    }
    else
    {
        ConjugateGradientMethodLinSolve(fluid, X, Xprev, b, A);
    }
}

void Fluids::advect(Fluid3D& fluid, std::vector<double>& D, const std::vector<double>& Dprev, const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& Z, const std::uint8_t b) const
{
    const std::uint64_t N = fluid.N; 
	const double dt = fluid.dt * N;

    std::uint64_t maxX = b == 1 ? fluid.N+1 : fluid.N+1;
    std::uint64_t maxY = b == 2 ? fluid.N+1 : fluid.N+1;

    /*
    std::cout << std::fixed << std::setprecision(2) << std::endl;
    std::cout << std::endl << "=== X ===" << std::endl;
    for (std::uint64_t j = 0; j < fluid.N+2; ++j)
    {
        for (std::uint64_t i = 0; i < fluid.N+2+1; ++i)
        {
            if (fluid.IX(i,j,0,1) < 10)
                std::cout << " ";
            if (fluid.IX(i,j,0,1) < 100)
                std::cout << " ";
            std::cout << fluid.IX(i,j,0,1) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << "=== Y ===" << std::endl;
    for (std::uint64_t j = 0; j < fluid.N+2+1; ++j)
    {
        for (std::uint64_t i = 0; i < fluid.N+2; ++i)
        {
            if (fluid.IX(i,j,0,2) < 10)
                std::cout << " ";
            if (fluid.IX(i,j,0,2) < 100)
                std::cout << " ";
            std::cout << fluid.IX(i,j,0,2) << " ";
        }
        std::cout << std::endl;
    }
    */

    for (std::uint64_t j = 1; j < maxY; ++j)
    {
        for (std::uint64_t i = 1; i < maxX; ++i)
        {
            const double xvel = b == 0 ? 0.5*(X[fluid.IX(i,j,0,1)]+X[fluid.IX(i+1,j,0,1)]) : X[fluid.IX(i,j,0,1)];
            const double yvel = b == 0 ? 0.5*(Y[fluid.IX(i,j,0,2)]+Y[fluid.IX(i,j+1,0,2)]) : Y[fluid.IX(i,j,0,2)];
            double x = std::clamp((i-dt*xvel) - (b == 2 ? 0 : 0.0), 0.5, N + 0.5);
            double y = std::clamp((j-dt*yvel) - (b == 1 ? 0 : 0.0), 0.5, N + 0.5);

            const std::uint16_t i0 = static_cast<std::uint16_t>(x);
            const std::uint16_t i1 = i0 + 1;
            const std::uint16_t j0 = static_cast<std::uint16_t>(y);
            const std::uint16_t j1 = j0 + 1;

            const double s1 = x	- i0;
            const double s0 = 1.0 - s1;
            const double t1 = y	- j0;
            const double t0 = 1.0 - t1;

            if (b == 0 && x > 2 && x < 4 && y > 2 && y < 4)
            {
                std::cout << b+0 << std::endl;
                std::cout << x << " " << y << std::endl;
                std::cout << i0 << " " << j0 << " " << fluid.IX(i0,j0,0,b) << std::endl;
                std::cout << i0 << " " << j1 << " " << fluid.IX(i0,j1,0,b) << std::endl;
                std::cout << i1 << " " << j0 << " " << fluid.IX(i1,j0,0,b) << std::endl;
                std::cout << i1 << " " << j1 << " " << fluid.IX(i1,j1,0,b) << std::endl;
                std::cout << std::endl;
            }

            double vA = Dprev[fluid.IX(i0,j0,0,b)];
            double vB = Dprev[fluid.IX(i0,j1,0,b)];
            double vC = Dprev[fluid.IX(i1,j0,0,b)];
            double vD = Dprev[fluid.IX(i1,j1,0,b)];

            D[fluid.IX(i,j,0,b)] = s0*(t0*vA+t1*vB)+s1*(t0*vC+t1*vD);
        }
    }
	setBnd(fluid, D, b);
}

void Fluids::GaussSeidelRelaxationLinSolve(const Fluid3D& fluid, std::vector<double>& X, const std::vector<double>& Xprev, const double a, const double c, std::uint8_t b) const
{
    std::uint64_t maxX = b == 1 ? fluid.N+2 : fluid.N+1;
    std::uint64_t maxY = b == 2 ? fluid.N+2 : fluid.N+1;

	const double cinv = 1.0/c;
	for (std::uint8_t l = 0; l < 64; ++l)
	{
        for (std::uint64_t j = 1; j < maxY; ++j)
        {
            for (std::uint64_t i = 1; i < maxX; ++i)
            {
                auto k = 0;
                X[fluid.IX(i,j,k,b)] =
                    (Xprev[fluid.IX(i,j,k,b)]
                        +a*(X[fluid.IX(i+1,j,k,b)]+X[fluid.IX(i-1,j,k,b)]+
                            X[fluid.IX(i,j+1,k,b)]+X[fluid.IX(i,j-1,k,b)]
                       ))*cinv;
            }
        }
		setBnd(fluid, X, b);
	}
}

void Fluids::applyPreconditioner(const std::uint64_t N, const Eigen::VectorXd& r, const Laplacian& A, Eigen::VectorXd& z, const Solver solver) const
{
    if (solver == CG)
    {
        z = r;
        return;
    }
    const std::uint64_t N2 = N*N; 
    const std::uint64_t N3 = N*N*N;

    // Solve Lq = r
    Eigen::VectorXd q = Eigen::VectorXd::Zero(z.size());
    for (std::int64_t n = 0; n < z.size(); ++n)
    {
        const std::uint64_t m = n % N2;
        const std::uint64_t i = m % N;
        const std::uint64_t j = m / N;
        const std::uint64_t k = n / N2;

        const std::uint64_t indmi = (i-1)+j*N+k*N2;
        const std::uint64_t indmj = i+(j-1)*N+k*N2;
        const std::uint64_t indmk = i+j*N+(k-1)*N2;

        const double a = i > 0 ? A.plusi.coeff(indmi) * A.precon.coeff(indmi) * q.coeff(indmi) : 0;
        const double b = j > 0 ? A.plusj.coeff(indmj) * A.precon.coeff(indmj) * q.coeff(indmj) : 0;
        const double c = k > 0 ? A.plusk.coeff(indmk) * A.precon.coeff(indmk) * q.coeff(indmk) : 0;

        const double t = r.coeff(n) - a - b - c;
        q.coeffRef(n) = t * A.precon.coeff(n);
    }

    // Solve L'z = q
    for (std::int64_t n = z.size()-1; n >= 0; --n)
    {
        const std::uint64_t m = n % N2;
        const std::uint64_t i = m % N;
        const std::uint64_t j = m / N;
        const std::uint64_t k = n / N2;

        const std::uint64_t indpi = (i+1)+j*N+k*N2;
        const std::uint64_t indpj = i+(j+1)*N+k*N2;
        const std::uint64_t indpk = i+j*N+(k+1)*N2;

        const double a = indpi < N3 ? z.coeff(indpi) : 0;
        const double b = indpj < N3 ? z.coeff(indpj) : 0;
        const double c = indpk < N3 ? z.coeff(indpk) : 0;

        const double prec = A.precon.coeff(n);
        const double t = q.coeff(n) - A.plusi.coeff(n) * prec * a
                                    - A.plusj.coeff(n) * prec * b
                                    - A.plusk.coeff(n) * prec * c;
        z.coeffRef(n) = t * prec;
    }
}

void Fluids::ConjugateGradientMethodLinSolve(const Fluid3D& fluid, std::vector<double>& X, const std::vector<double>& Xprev, const std::uint8_t bs, const Laplacian& A)
{
    const std::uint64_t N = fluid.N+1;
    const std::uint64_t N2 = (N+1)*(N+1);
    const std::uint32_t diagSize = A.diag.size();
    Eigen::VectorXd x(diagSize);
    Eigen::VectorXd b(diagSize);

    // Filling matrices
    #pragma omp parallel for schedule(static)
    for (std::uint64_t n = 0; n < diagSize; ++n)
    {
        const std::uint64_t m = n % N2;
        const std::uint64_t i = m % N;
        const std::uint64_t j = m / N;
        const std::uint64_t k = n / N2;
        b.coeffRef(n) = Xprev[fluid.IX(i+1,j+1,k+1,bs)];
    }

    // Solving Ap = b
    Eigen::VectorXd r = b;
    if (r.isZero(0))
    {
        X = Xprev;
        return;
    }
    Eigen::VectorXd p = Eigen::VectorXd::Zero(diagSize);
    Eigen::VectorXd z = p;
    applyPreconditioner(N, r, A, z, fluid.solver);
    Eigen::VectorXd s = z;
    double sig = z.dot(r);

    for (std::uint32_t i = 0; i < b.size(); ++i)
    {
        z = A.A * s;
        const double alpha = sig / s.dot(z);
        p = p + alpha * s;
        r = r - alpha * z;
        if (r.lpNorm<Eigen::Infinity>() < 10e-5)
        {
            break;
        }
        applyPreconditioner(N, r, A, z, fluid.solver);
        const double signew = z.dot(r);
        const double beta = signew / sig;
        s = z + beta * s;
        sig = signew;
    }

    // Write the results
    //#pragma omp parallel for schedule(static)
    for (std::uint64_t n = 0; n < diagSize; ++n)
    {
        const std::uint64_t m = n % N2;
        const std::uint64_t i = m % N;
        const std::uint64_t j = m / N;
        const std::uint64_t k = n / N2;
        X[fluid.IX(i+1,j+1,k+1,bs)] = p.coeff(n);
    }


    setBnd(fluid, X, bs);
}

void Fluids::project(const Fluid3D& fluid, std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Z, std::vector<double>& p, std::vector<double>& div)
{
    const double h = 1.0/fluid.N;

    for (std::uint16_t j = 1; j < fluid.N+1; ++j)
    {
        for (std::uint16_t i = 1; i < fluid.N+1; ++i)
        {
            auto k = 0;
            div[fluid.IX(i,j,k,0)] = - 0.5 * (
                        (X[fluid.IX(i+1,j,k,1)]-X[fluid.IX(i,j,k,1)])+
                        (Y[fluid.IX(i,j+1,k,2)]-Y[fluid.IX(i,j,k,2)])
                        )*h;
        }
    }
    for (std::uint16_t j = 0; j < fluid.N+2; ++j)
    {
        for (std::uint16_t i = 0; i < fluid.N+2; ++i)
        {
            auto k = 0;
            p[fluid.IX(i,j,k,0)] = 0;
        }
    }
    //p.clear();

	setBnd(fluid, div, 0);

    if (fluid.solver == GAUSS_SEIDEL)
    {
        GaussSeidelRelaxationLinSolve(fluid, p, div, 1, fluid.is2D ? 4 : 6, 0);
    }
    else
    {
        ConjugateGradientMethodLinSolve(fluid, p, div, 0, fluid.laplacianProject);
    }

    for (std::uint16_t j = 1; j < fluid.N+1; ++j)
    {
        for (std::uint16_t i = 1; i < fluid.N+1; ++i)
        {
            auto k = 0;
            X[fluid.IX(i,j,k,1)] -= 0.5*fluid.N*(p[fluid.IX(i,j,k,0)]-p[fluid.IX(i-1,j,k,0)]);
            Y[fluid.IX(i,j,k,2)] -= 0.5*fluid.N*(p[fluid.IX(i,j,k,0)]-p[fluid.IX(i,j-1,k,0)]);
        }
    }

	setBnd(fluid, X, 1);
	setBnd(fluid, Y, 2);
}

void Fluids::setBnd(const Fluid3D& fluid, std::vector<double>& X, const std::uint8_t b) const
{
    std::uint64_t maxX = b == 1 ? fluid.N+2 : fluid.N+1;
    std::uint64_t maxY = b == 2 ? fluid.N+2 : fluid.N+1;

    for (std::uint16_t i = 1; i <= fluid.N+1; ++i)
    {
        if (b != 1 || (b == 1 && i <= fluid.N))
        {
            X[fluid.IX(0,i,0,b)] = b == 1 ? -X[fluid.IX(1,i,0,b)] : X[fluid.IX(1,i,0,b)];
            X[fluid.IX(maxX,i,0,b)] = b == 1 ? -X[fluid.IX(maxX-1,i,0,b)] : X[fluid.IX(maxX-1,i,0,b)];
        }
        if (b != 2 || (b == 2 && i <= fluid.N))
        {
            X[fluid.IX(i,0,0,b)] = b == 2 ? -X[fluid.IX(i,1,0,b)] : X[fluid.IX(i,1,0,b)];
            X[fluid.IX(i,maxY,0,b)] = b == 2 ? -X[fluid.IX(i,maxY-1,0,b)] : X[fluid.IX(i,maxY-1,0,b)];
        }
    }

    X[fluid.IX(0,0,0,b)]        = 0.5*(X[fluid.IX(1,0,0,b)]+X[fluid.IX(0,1,0,b)]);
    X[fluid.IX(0,maxY,0,b)]     = 0.5*(X[fluid.IX(1,maxY,0,b)]+X[fluid.IX(0,maxY-1,0,b)]);
    X[fluid.IX(maxX,0,0,b)]     = 0.5*(X[fluid.IX(maxX-1,0,0,b)]+X[fluid.IX(maxX,1,0,b)]);
    X[fluid.IX(maxX,maxY,0,b)]  = 0.5*(X[fluid.IX(maxX,maxY-1,0,b)]+X[fluid.IX(maxX-1,maxY,0,b)]);
}

void Fluids::updateRender(Fluid3D& fluid)
{
    const std::uint64_t N = fluid.N;
    const std::uint64_t N32 = fluid.is2D ? (N)*(N) : (N+2)*(N+2)*(N+2);
	std::vector<std::uint8_t> texture(fluid.is2D ? N32*3 : N32, 0);

    std::uint64_t it = 0;
    for (std::uint16_t j = 1; j < fluid.N+1; ++j)
    {
        for (std::uint16_t i = 1; i < fluid.N+1; ++i)
        {
		    texture[it*3] = static_cast<std::uint8_t>(std::clamp(fluid.substanceField[fluid.IX(i,j,0,0)], 0.0, 255.0));
		    texture[it*3+1] = static_cast<std::uint8_t>(std::clamp(fluid.substanceField[fluid.IX(i,j,0,0)], 0.0, 255.0));
		    texture[it*3+2] = static_cast<std::uint8_t>(std::clamp(fluid.substanceField[fluid.IX(i,j,0,0)], 0.0, 255.0));
            it++;
            /*
            double implicit = fluid.implicitFunctionField[i];
            texture[i*3] = 0;
            texture[i*3+1] = 0;
            texture[i*3+2] = 0;
            std::uint64_t p = 50;
            if (implicit < 0)
            {
                texture[i*3+2] = std::clamp(-implicit*p, 0.0, 255.0);
            }
            else
            {
                texture[i*3] = std::clamp(implicit*p, 0.0, 255.0);
            }
            */
            /*
            const std::uint64_t m = i % N32;
            const std::uint64_t x = m % (N + 2);
            const std::uint64_t y = m / (N + 2);
            double gO = gradLength(fluid, fluid.implicitFunctionFieldPrev, x, y);
            texture[i*3] = std::clamp(gO*128, 0.0, 255.0);
            texture[i*3+1] = 0;
            texture[i*3+2] = 0;
            */
        }
	}

	const auto& textureGL = gCoordinator.GetComponent<Material>(fluid.entity).texture;
    if (!fluid.is2D)
    {
	    _renderer->initTexture3D(texture, textureGL);
    }
    else
    {
	    _renderer->initTexture2D(texture, textureGL);
    }
}

void Fluids::writeVolumeFile(Fluid3D& fluid, std::uint64_t iteration)
{
    const std::uint64_t N = fluid.N;
    const std::uint64_t N2 = N*N;
    const std::uint64_t N22 = (N+2)*(N+2);
    const std::uint64_t N3 = N*N*N;
    std::string path = "result/";
    path += std::to_string(iteration);
    path += ".vol";
    std::ofstream file (path, std::ios::binary | std::ofstream::trunc);
    file << 'V';
    file << 'O';
    file << 'L';
    write(file, (uint8_t)3); // Version
    write(file, (int32_t)1); // Type
    std::uint16_t n = fluid.N;
    write(file, (int32_t)n);
    write(file, (int32_t)n);
    write(file, (int32_t)n);
    write(file, (int32_t)1); // Nb channels
    float xmin = -0.5f;
    float ymin = -0.5f;
    float zmin = -0.5f;
    float xmax = 0.5f;
    float ymax = 0.5f;
    float zmax = 0.5f;
    write(file, xmin);
    write(file, ymin);
    write(file, zmin);
    write(file, xmax);
    write(file, ymax);
    write(file, zmax);

    for (std::uint64_t n = 0; n < N3; ++n)
    {
        const std::uint32_t m = n % N2;
        const std::uint16_t i = m % N + 1;
        const std::uint16_t j = m / N + 1;
        const std::uint16_t k = n / N2 + 1;
        const std::uint64_t ind = i+j*(N+2)+k*N22;
        float value = std::clamp(fluid.substanceField[ind], 0.0, 255.0);
        write(file, value);
    }

    file.close();
}
