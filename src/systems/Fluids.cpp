#include "Fluids.h"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <Eigen/src/SparseLU/SparseLU.h>
#include <cmath>
#include <cstdint>
#include <omp.h>

void Fluids::init(std::shared_ptr<Renderer> renderer)
{
	_renderer = renderer;
    std::srand(std::time(nullptr));
    omp_set_num_threads(0);
}

void Fluids::update(std::uint32_t iteration)
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

        float p = 64;
        float z = 512;

        if (iteration >= 0)
        {
            for (std::uint32_t i = 0; i < 16; ++i)
            {
                for (std::uint32_t j = 0; j < 16; ++j)
                {
                    fluid.velocityFieldZ[fluid.IX(N+i-(i/2), N+j-(j/2), N)] = -z;
                    fluid.substanceField[fluid.IX(N+i-(i/2), N+j-(j/2), N)] = p;
                    /*
                    fluid.velocityFieldY[fluid.IX(N+i-(i/2), fluid.N-3, N+j-(j/2))] = -z;
                    fluid.substanceField[fluid.IX(N+i-(i/2), fluid.N-3, N+j-(j/2))] = p;
                    */

                }
            }
        }

        /*
        float e = 1024*9;
        float u = 1024*9;

	    fluid.velocityFieldX[fluid.IX(N, fluid.N-2, 2)] = -e;
	    fluid.velocityFieldZ[fluid.IX(N, fluid.N-2, 2)] = e*1.5;
	    fluid.velocityFieldX[fluid.IX(N, fluid.N-2, fluid.N-2)] = e;
	    fluid.velocityFieldZ[fluid.IX(N, fluid.N-2, fluid.N-2)] = -e*1.5;
	    fluid.velocityFieldZ[fluid.IX(2, fluid.N-2, N)] = e;
	    fluid.velocityFieldX[fluid.IX(2, fluid.N-2, N)] = e*1.5;
	    fluid.velocityFieldZ[fluid.IX(fluid.N-2, fluid.N-2, N)] = -e;
	    fluid.velocityFieldX[fluid.IX(fluid.N-2, fluid.N-2, N)] = -e*1.5;


	    fluid.velocityFieldX[fluid.IX(N, 2, 2)] = -e;
	    fluid.velocityFieldX[fluid.IX(N, 2, fluid.N-2)] = e;
	    fluid.velocityFieldZ[fluid.IX(2, 2, N)] = e;
	    fluid.velocityFieldZ[fluid.IX(fluid.N-2, 2, N)] = -e;
	    fluid.velocityFieldY[fluid.IX(N, 2, 2)] = u;
	    fluid.velocityFieldY[fluid.IX(N, 2, fluid.N-2)] = u;
	    fluid.velocityFieldY[fluid.IX(2, 2, N)] = u;
	    fluid.velocityFieldY[fluid.IX(fluid.N-2, 2, N)] = u;

	    fluid.velocityFieldX[fluid.IX(N, N, 2)] = -e;
	    fluid.velocityFieldX[fluid.IX(N, N, fluid.N-2)] = e;
	    fluid.velocityFieldZ[fluid.IX(2, N, N)] = e;
	    fluid.velocityFieldZ[fluid.IX(fluid.N-2, N, N)] = -e;
	    fluid.velocityFieldY[fluid.IX(N, N, 2)] = u;
	    fluid.velocityFieldY[fluid.IX(N, N, fluid.N-2)] = u;
	    fluid.velocityFieldY[fluid.IX(2, N, N)] = u;
	    fluid.velocityFieldY[fluid.IX(fluid.N-2, N, N)] = u;
        */

        /* End Testings */

		Vstep(fluid);
		Sstep(fluid);
		updateRender(fluid);
	}
}


void Fluids::Vstep(Fluid3D& fluid)
{
    /*
	addSource(fluid, fluid.velocityFieldX, fluid.velocityFieldPrevX);
	addSource(fluid, fluid.velocityFieldY, fluid.velocityFieldPrevY);
	addSource(fluid, fluid.velocityFieldZ, fluid.velocityFieldPrevZ);
    */

    /*
    float a = fluid.dt * fluid.viscosity * fluid.N;
    GaussSeidelRelaxationLinSolve(fluid, fluid.velocityFieldPrevX, fluid.velocityFieldX, a, 1+6*a, 1);
    GaussSeidelRelaxationLinSolve(fluid, fluid.velocityFieldPrevY, fluid.velocityFieldY, a, 1+6*a, 2);
    GaussSeidelRelaxationLinSolve(fluid, fluid.velocityFieldPrevZ, fluid.velocityFieldZ, a, 1+6*a, 3);
    */
	diffuse(fluid, fluid.velocityFieldPrevX, fluid.velocityFieldX, 1, fluid.laplacianViscosity);
	diffuse(fluid, fluid.velocityFieldPrevY, fluid.velocityFieldY, 2, fluid.laplacianViscosity);
	diffuse(fluid, fluid.velocityFieldPrevZ, fluid.velocityFieldZ, 3, fluid.laplacianViscosity);

	project(fluid, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, fluid.velocityFieldX, fluid.velocityFieldY);

	advect(fluid, fluid.velocityFieldX, fluid.velocityFieldPrevX, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, 1);
	advect(fluid, fluid.velocityFieldY, fluid.velocityFieldPrevY, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, 2);
	advect(fluid, fluid.velocityFieldZ, fluid.velocityFieldPrevZ, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY, fluid.velocityFieldPrevZ, 3);

	project(fluid, fluid.velocityFieldX, fluid.velocityFieldY, fluid.velocityFieldZ, fluid.velocityFieldPrevX, fluid.velocityFieldPrevY);
}

void Fluids::Sstep(Fluid3D& fluid)
{
	//addSource(fluid, fluid.substanceField, fluid.substanceFieldPrev);
    /*
    float a = fluid.dt * fluid.diffusion * fluid.N;
    GaussSeidelRelaxationLinSolve(fluid, fluid.substanceFieldPrev, fluid.substanceField, a, 1+6*a, 0);
    */
	diffuse(fluid, fluid.substanceFieldPrev, fluid.substanceField, 0, fluid.laplacianDiffuse);
	advect(fluid, fluid.substanceField, fluid.substanceFieldPrev, fluid.velocityFieldX, fluid.velocityFieldY, fluid.velocityFieldZ, 0);
}

void Fluids::addSource(const Fluid3D& fluid, std::vector<double>& X, const std::vector<double>& S) const
{
    #pragma omp for
	for (std::uint32_t i = 0; i < (fluid.N+2)*(fluid.N+2)*(fluid.N+2); ++i)
	{
		X[i] += fluid.dt * S[i];
	}
}

void Fluids::diffuse(const Fluid3D& fluid, std::vector<double>& X, const std::vector<double>& Xprev, const std::uint8_t b, const Laplacian& A)
{
    ConjugateGradientMethodLinSolve(fluid, X, Xprev, b, A);
}

void Fluids::advect(Fluid3D& fluid, std::vector<double>& D, const std::vector<double>& Dprev, const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& Z, const std::uint8_t b) const
{
    const std::uint64_t N = fluid.N; 
    const std::uint64_t N2 = N*N; 
    const std::uint64_t N22 = (N+2)*(N+2); 
    const std::uint64_t N3 = N*N*N; 
	const double dt = fluid.dt * N;

    #pragma omp parallel for
    for (std::uint64_t n = 0; n < N3; ++n)
    {
        const std::uint64_t m = n % N2;
        const std::uint64_t i = m % N + 1;
        const std::uint64_t j = m / N + 1;
        const std::uint64_t k = n / N2 + 1;
        const std::uint64_t ind = i+j*(N+2)+k*N22;

        const double x = std::clamp(i-dt*X[ind], 0.5, N + 0.5);
        const double y = std::clamp(j-dt*Y[ind], 0.5, N + 0.5);
        const double z = std::clamp(k-dt*Z[ind], 0.5, N + 0.5);
        const std::uint16_t i0 = static_cast<std::uint16_t>(x);
        const std::uint16_t i1 = i0 + 1;
        const std::uint16_t j0 = static_cast<std::uint16_t>(y);
        const std::uint16_t j1 = j0 + 1;
        const std::uint16_t k0 = static_cast<std::uint16_t>(z);
        const std::uint16_t k1 = k0 + 1;

        const double s1 = x	- i0;
        const double s0 = 1.0 - s1;
        const double t1 = y	- j0;
        const double t0 = 1.0 - t1;
        const double u1 = z	- k0;
        const double u0 = 1.0 - u1;

        D[ind] =
            s0*(t0*(u0*Dprev[fluid.IX(i0,j0,k0)]
                    +u1*Dprev[fluid.IX(i0,j0,k1)])
                +(t1*(u0*Dprev[fluid.IX(i0,j1,k0)]
                    +u1*Dprev[fluid.IX(i0,j1,k1)])))
            +s1*(t0*(u0*Dprev[fluid.IX(i1,j0,k0)]
                    +u1*Dprev[fluid.IX(i1,j0,k1)])
                +(t1*(u0*Dprev[fluid.IX(i1,j1,k0)]
                    +u1*Dprev[fluid.IX(i1,j1,k1)])));
    }

    // Reverse advection to calculate errors made, than correct the first advection to reduce the errors
    #pragma omp parallel for
    for (std::uint64_t n = 0; n < N3; ++n)
    {
        const std::uint64_t m = n % N2;
        const std::uint64_t i = m % N + 1;
        const std::uint64_t j = m / N + 1;
        const std::uint64_t k = n / N2 + 1;
        const std::uint64_t ind = i+j*(N+2)+k*N22;

        // Backward to get neighbours for clamping
        double x = std::clamp(i-dt*X[ind], 0.5, N + 0.5);
        double y = std::clamp(j-dt*Y[ind], 0.5, N + 0.5);
        double z = std::clamp(k-dt*Z[ind], 0.5, N + 0.5);
        std::uint16_t i0 = static_cast<std::uint16_t>(x);
        std::uint16_t i1 = i0 + 1;
        std::uint16_t j0 = static_cast<std::uint16_t>(y);
        std::uint16_t j1 = j0 + 1;
        std::uint16_t k0 = static_cast<std::uint16_t>(z);
        std::uint16_t k1 = k0 + 1;

        const double top = std::max({D[fluid.IX(i0,j0,k0)],D[fluid.IX(i0,j0,k1)],D[fluid.IX(i0,j1,k0)],D[fluid.IX(i0,j1,k1)],D[fluid.IX(i1,j0,k0)],D[fluid.IX(i1,j0,k1)],D[fluid.IX(i1,j1,k0)],D[fluid.IX(i1,j1,k1)]});
        const double bot = std::min({D[fluid.IX(i0,j0,k0)],D[fluid.IX(i0,j0,k1)],D[fluid.IX(i0,j1,k0)],D[fluid.IX(i0,j1,k1)],D[fluid.IX(i1,j0,k0)],D[fluid.IX(i1,j0,k1)],D[fluid.IX(i1,j1,k0)],D[fluid.IX(i1,j1,k1)]});

        // Forward after backward to get error
        x = std::clamp(i+dt*X[ind], 0.5, N + 0.5);
        y = std::clamp(j+dt*Y[ind], 0.5, N + 0.5);
        z = std::clamp(k+dt*Z[ind], 0.5, N + 0.5);
        i0 = static_cast<std::uint16_t>(x);
        i1 = i0 + 1;
        j0 = static_cast<std::uint16_t>(y);
        j1 = j0 + 1;
        k0 = static_cast<std::uint16_t>(z);
        k1 = k0 + 1;

        const double s1 = x	- i0;
        const double s0 = 1.0 - s1;
        const double t1 = y	- j0;
        const double t0 = 1.0 - t1;
        const double u1 = z	- k0;
        const double u0 = 1.0 - u1;

        const double back = 
            s0*(t0*(u0*D[fluid.IX(i0,j0,k0)]
                    +u1*D[fluid.IX(i0,j0,k1)])
                +(t1*(u0*D[fluid.IX(i0,j1,k0)]
                    +u1*D[fluid.IX(i0,j1,k1)])))
            +s1*(t0*(u0*D[fluid.IX(i1,j0,k0)]
                    +u1*D[fluid.IX(i1,j0,k1)])
                +(t1*(u0*D[fluid.IX(i1,j1,k0)]
                    +u1*D[fluid.IX(i1,j1,k1)])));
        D[ind] = std::clamp(D[ind] + 0.5 * (Dprev[ind] - back), bot, top);
	}
	setBnd(fluid, D, b);
}

void Fluids::GaussSeidelRelaxationLinSolve(const Fluid3D& fluid, std::vector<double>& X, std::vector<double>& Xprev, float a, float c, std::uint8_t b) const
{
    exit(0);
	float cinv = 1.0f/c;
	for (std::uint32_t l = 0; l < 8; ++l)
	{
		for (std::uint32_t k = 1; k <= fluid.N; ++k)
        {
			for (std::uint32_t j = 1; j <= fluid.N; ++j)
            {
				for (std::uint32_t i = 1; i <= fluid.N; ++i)
                {
					X[fluid.IX(i,j,k)] =
                        (Xprev[fluid.IX(i,j,k)]
							+a*(X[fluid.IX(i+1,j,k)]+
								X[fluid.IX(i-1,j,k)]+
								X[fluid.IX(i,j+1,k)]+
								X[fluid.IX(i,j-1,k)]+
								X[fluid.IX(i,j,k+1)]+
								X[fluid.IX(i,j,k-1)]
						   ))*cinv;
				}
			}
		}
		setBnd(fluid, X, b);
	}
}

void Fluids::applyPreconditioner(const std::uint64_t N, const Eigen::VectorXd& r, const Laplacian& A, Eigen::VectorXd& z) const
{
    const std::uint64_t N2 = N*N; 
    const std::uint64_t N3 = N*N*N; 

    // Solve Lq = r
    Eigen::VectorXd q = Eigen::VectorXd::Zero(z.size());
    for (std::uint32_t k = 0; k < N; ++k)
    {
		for (std::uint32_t j = 0; j < N; ++j)
        {
			for (std::uint32_t i = 0; i < N; ++i)
            {
                std::uint32_t ind = i+j*N+k*N2;
                if (ind < z.size())
                {
                    std::uint32_t indmi = (i-1)+j*N+k*N2;
                    std::uint32_t indmj = i+(j-1)*N+k*N2;
                    std::uint32_t indmk = i+j*N+(k-1)*N2;

                    double a = i > 0 ? A.plusi.coeff(indmi) * A.precon.coeff(indmi) * q.coeff(indmi) : 0;
                    double b = j > 0 ? A.plusj.coeff(indmj) * A.precon.coeff(indmj) * q.coeff(indmj) : 0;
                    double c = k > 0 ? A.plusk.coeff(indmk) * A.precon.coeff(indmk) * q.coeff(indmk) : 0;

                    double t = r.coeff(ind) - a - b - c;
                    q.coeffRef(ind) = t * A.precon.coeff(ind);
                }
            }
        }
    }

    // Solve L.transpose()z = q
    for (std::int32_t k = N-1; k >= 0; --k)
    {
		for (std::int32_t j = N-1; j >= 0; --j)
        {
			for (std::int32_t i = N-1; i >= 0; --i)
            {
                std::uint32_t ind   = i+j*N+k*N2;
                if (ind < z.size())
                {
                    std::uint32_t indpi = (i+1)+j*N+k*N2;
                    std::uint32_t indpj = i+(j+1)*N+k*N2;
                    std::uint32_t indpk = i+j*N+(k+1)*N2;

                    double a = indpi < N3 ? z.coeff(indpi) : 0;
                    double b = indpj < N3 ? z.coeff(indpj) : 0;
                    double c = indpk < N3 ? z.coeff(indpk) : 0;

                    double prec = A.precon.coeff(ind);
                    double t = q[ind]   - A.plusi.coeff(ind) * prec * a
                                        - A.plusj.coeff(ind) * prec * b
                                        - A.plusk.coeff(ind) * prec * c;
                    z.coeffRef(ind) = t * prec;
                }
            }
        }
    }
}

void Fluids::ConjugateGradientMethodLinSolve(const Fluid3D& fluid, std::vector<double>& X, const std::vector<double>& Xprev, const std::uint8_t bs, const Laplacian& A)
{
    const std::uint64_t N = fluid.N;
    const std::uint64_t N2 = N*N;
    const std::uint32_t diagSize = A.diag.size();
    Eigen::VectorXd x(diagSize);
    Eigen::VectorXd b(diagSize);

    // Filling matrices
    #pragma omp parallel for
    for (std::uint64_t n = 0; n < diagSize; ++n)
    {
        const std::uint64_t m = n % N2;
        const std::uint64_t i = m % N;
        const std::uint64_t j = m / N;
        const std::uint64_t k = n / N2;
        b.coeffRef(n) = Xprev[fluid.IX(i+1,j+1,k+1)];
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
    applyPreconditioner(fluid.N, r, A, z);
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
        applyPreconditioner(fluid.N, r, A, z);
        const double signew = z.dot(r);
        const double beta = signew / sig;
        s = z + beta * s;
        sig = signew;
    }

    // Write the results
    #pragma omp parallel for
    for (std::uint64_t n = 0; n < diagSize; ++n)
    {
        const std::uint64_t m = n % N2;
        const std::uint64_t i = m % N;
        const std::uint64_t j = m / N;
        const std::uint64_t k = n / N2;
        X[fluid.IX(i+1,j+1,k+1)] = p.coeff(n);
    }

    setBnd(fluid, X, bs);
}

void Fluids::project(const Fluid3D& fluid, std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Z, std::vector<double>& p, std::vector<double>& div)
{
    const double h = 1.0/fluid.N;
    const std::uint64_t N = fluid.N;
    const std::uint64_t N2 = N*N;
    const std::uint64_t N22 = (N+2)*(N+2);
    const std::uint64_t N3 = N*N*N;

    #pragma omp parallel for
    for (std::uint64_t n = 0; n < N3; ++n)
    {
        const std::uint64_t m = n % N2;
        const std::uint64_t i = m % N + 1;
        const std::uint64_t j = m / N + 1;
        const std::uint64_t k = n / N2 + 1;
        const std::uint64_t ind = i+j*(N+2)+k*N22;

        div[ind] = -(1.0/3.0)*
                    ((X[fluid.IX(i+1,j,k)]-X[fluid.IX(i-1,j,k)])*h+
                    (Y[fluid.IX(i,j+1,k)]-Y[fluid.IX(i,j-1,k)])*h+
                    (Z[fluid.IX(i,j,k+1)]-Z[fluid.IX(i,j,k-1)])*h);
    }
    p.clear();

	setBnd(fluid, div, 0);
	setBnd(fluid, p, 0);

    ConjugateGradientMethodLinSolve(fluid, p, div, 0, fluid.laplacianProject);
    //GaussSeidelRelaxationLinSolve(fluid, p, div, 1, 6, 0);

    #pragma omp parallel for
    for (std::uint64_t n = 0; n < N3; ++n)
    {
        const std::uint32_t m = n % N2;
        const std::uint16_t i = m % N + 1;
        const std::uint16_t j = m / N + 1;
        const std::uint16_t k = n / N2 + 1;
        const std::uint64_t ind = i+j*(N+2)+k*N22;

        X[ind] -= 0.5*N*(p[fluid.IX(i+1,j,k)]-p[fluid.IX(i-1,j,k)]);
        Y[ind] -= 0.5*N*(p[fluid.IX(i,j+1,k)]-p[fluid.IX(i,j-1,k)]);
        Z[ind] -= 0.5*N*(p[fluid.IX(i,j,k+1)]-p[fluid.IX(i,j,k-1)]);
    }

	setBnd(fluid, X, 1);
	setBnd(fluid, Y, 2);
	setBnd(fluid, Z, 3);
}

void Fluids::setBnd(const Fluid3D& fluid, std::vector<double>& X, const std::uint8_t b) const
{
    const std::uint64_t N = fluid.N;
    const std::uint64_t N2 = N*N;
    for (std::uint32_t n = 0; n < N2; ++n)
    {
        const std::uint32_t m = n % N2;
        const std::uint16_t i = m % N + 1;
        const std::uint16_t j = m / N + 1;

        X[fluid.IX(i,j,0)]	= b == 3 ? -X[fluid.IX(i,j,1)]		: X[fluid.IX(i,j,1)];
        X[fluid.IX(i,j,N)]  = b == 3 ? -X[fluid.IX(i,j,N-1)]	: X[fluid.IX(i,j,N-1)];

        X[fluid.IX(i,0,j)]	= b == 2 ? -X[fluid.IX(i,1,j)]		: X[fluid.IX(i,1,j)];
        X[fluid.IX(i,N,j)]  = b == 2 ? -X[fluid.IX(i,N-1,j)]	: X[fluid.IX(i,N-1,j)];

        X[fluid.IX(0,i,j)]	= b == 1 ? -X[fluid.IX(1,i,j)]	    : X[fluid.IX(1,i,j)];
        X[fluid.IX(N,i,j)]  = b == 1 ? -X[fluid.IX(N-1,i,j)]	: X[fluid.IX(N-1,i,j)];
    }

    // Edges
    for (std::uint32_t i = 1; i <= N; ++i)
    {
        X[fluid.IX(i,0,0)]  = 0.5*(X[fluid.IX(i,1,0)]   +X[fluid.IX(i,0,1)]);
        X[fluid.IX(i,N,0)]  = 0.5*(X[fluid.IX(i,N-1,0)] +X[fluid.IX(i,N,1)]);
        X[fluid.IX(i,0,N)]  = 0.5*(X[fluid.IX(i,0,N-1)] +X[fluid.IX(i,1,N)]);
        X[fluid.IX(i,N,N)]  = 0.5*(X[fluid.IX(i,N-1,N)] +X[fluid.IX(i,N,N-1)]);

        X[fluid.IX(0,i,0)]  = 0.5*(X[fluid.IX(1,i,0)]   +X[fluid.IX(0,i,1)]);
        X[fluid.IX(N,i,0)]  = 0.5*(X[fluid.IX(N-1,i,0)] +X[fluid.IX(N,i,1)]);
        X[fluid.IX(0,i,N)]  = 0.5*(X[fluid.IX(0,i,N-1)] +X[fluid.IX(1,i,N)]);
        X[fluid.IX(N,i,N)]  = 0.5*(X[fluid.IX(N-1,i,N)] +X[fluid.IX(N,i,0)]);

        X[fluid.IX(0,0,i)]  = 0.5*(X[fluid.IX(0,1,i)]   +X[fluid.IX(1,0,i)]);
        X[fluid.IX(0,N,i)]  = 0.5*(X[fluid.IX(0,N-1,i)] +X[fluid.IX(1,N,i)]);
        X[fluid.IX(N,0,i)]  = 0.5*(X[fluid.IX(N-1,0,i)] +X[fluid.IX(N,1,i)]);
        X[fluid.IX(N,N,i)]  = 0.5*(X[fluid.IX(N,N-1,i)] +X[fluid.IX(N-1,N,i)]);
    }

	// Corners
	X[fluid.IX(0,0,0)]	= (1.0/3.0)*(X[fluid.IX(1,0,0)] +X[fluid.IX(0,1,0)]     +X[fluid.IX(0,0,1)]);
	X[fluid.IX(0,N,0)]	= (1.0/3.0)*(X[fluid.IX(1,N,0)] +X[fluid.IX(0,N-1,0)]   +X[fluid.IX(0,N,1)]);
	X[fluid.IX(0,0,N)]  = (1.0/3.0)*(X[fluid.IX(1,0,N)] +X[fluid.IX(0,1,N)]     +X[fluid.IX(0,0,N+1)]);
	X[fluid.IX(0,N,N)]  = (1.0/3.0)*(X[fluid.IX(1,N,N)] +X[fluid.IX(0,N-1,N)]   +X[fluid.IX(0,N,N-1)]);

	X[fluid.IX(N,0,0)]	= (1.0/3.0)*(X[fluid.IX(N-1,0,0)]   +X[fluid.IX(N,1,0)]     +X[fluid.IX(N,0,1)]);
	X[fluid.IX(N,N,0)]	= (1.0/3.0)*(X[fluid.IX(N-1,N,0)]   +X[fluid.IX(N,N-1,0)]   +X[fluid.IX(N,N,1)]);
	X[fluid.IX(N,0,N)]	= (1.0/3.0)*(X[fluid.IX(N-1,0,N)]   +X[fluid.IX(N,1,N)]     +X[fluid.IX(N,0,N-1)]);
	X[fluid.IX(N,N,N)]  = (1.0/3.0)*(X[fluid.IX(N-1,N,N)]   +X[fluid.IX(N,N-1,N)]   +X[fluid.IX(N,N,N-1)]);
}

void Fluids::updateRender(Fluid3D& fluid)
{
    const std::uint64_t N = fluid.N;
    const std::uint64_t N32 = (N+2)*(N+2)*(N+2);
	std::vector<std::uint8_t> texture(N32, 0);

    #pragma omp parallel for
	for (std::uint64_t i = 0; i < N32; ++i)
	{
		texture[i] = static_cast<std::uint8_t>(std::clamp(fluid.substanceField[i], 0.0, 255.0));
	}

	const auto& textureGL = gCoordinator.GetComponent<Material>(fluid.entity).texture;
	_renderer->initTexture3D(texture, textureGL);
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
