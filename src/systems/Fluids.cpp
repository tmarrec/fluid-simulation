#include "Fluids.h"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <Eigen/src/SparseLU/SparseLU.h>
#include <cmath>
#include <cstdint>
#include <omp.h>

#ifdef DEBUG_GUI
float debugViscosity = 0;
float debugDiffusion = 0;
float debugDt = 0;
float debugAbsorption = 0;
float debugLightIntensity[3] = {0, 0, 0};
int debugN = 0;
#endif

void Fluids::init(std::shared_ptr<Renderer> renderer)
{
	_renderer = renderer;
    std::srand(std::time(nullptr));
}

void Fluids::reset([[maybe_unused]] bool force)
{
#ifdef DEBUG_GUI
	for (auto const& entity : mEntities)
	{
		auto& fluid = gCoordinator.GetComponent<Fluid3D>(entity);
        bool rebuild = debugDt != fluid.dt || debugViscosity != fluid.viscosity || debugN != (int)fluid.N || fluid.diffusion != debugDiffusion;
        if (force || rebuild)
        {
            fluid = Fluid3D();
            fluid.entity = entity;
            fluid.N = debugN;
            fluid.viscosity = debugViscosity;
            fluid.diffusion = debugDiffusion;
            fluid.dt = debugDt;
            fluid.init();
            _initCG = false;
        }
		auto& material = gCoordinator.GetComponent<Material>(entity);
        material.absorption = debugAbsorption;
        material.lightIntensity.x = debugLightIntensity[0];
        material.lightIntensity.y = debugLightIntensity[1];
        material.lightIntensity.z = debugLightIntensity[2];
    }
#endif
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
            for (std::uint32_t i = 0; i < 12; ++i)
            {
                for (std::uint32_t j = 0; j < 12; ++j)
                {
                    fluid.velocityFieldY[fluid.IX(N+i-(i/2), 3, N+j-(j/2))] = z;
                    if (iteration >= 0)
                    {
                        fluid.substanceField[fluid.IX(N+i-(i/2), 3, N+j-(j/2))] = p;
                    }
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
	double dt0 = fluid.dt * fluid.N;
    // First backward advection
	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t j = 1; j <= fluid.N; ++j)
		{
			for (std::uint32_t i = 1; i <= fluid.N; ++i)
			{
				double x = i-dt0*X[fluid.IX(i,j,k)];
				double y = j-dt0*Y[fluid.IX(i,j,k)];
				double z = k-dt0*Z[fluid.IX(i,j,k)];
				x = std::clamp(x, 0.5, fluid.N + 0.5);
				y = std::clamp(y, 0.5, fluid.N + 0.5);
				z = std::clamp(z, 0.5, fluid.N + 0.5);

				x = std::clamp(x, 0.5, fluid.N + 0.5);
				std::uint32_t i0 = static_cast<std::uint32_t>(x);
				std::uint32_t i1 = i0 + 1;
				y = std::clamp(y, 0.5, fluid.N + 0.5);
				std::uint32_t j0 = static_cast<std::uint32_t>(y);
				std::uint32_t j1 = j0 + 1;
				z = std::clamp(z, 0.5, fluid.N + 0.5);
				std::uint32_t k0 = static_cast<std::uint32_t>(z);
				std::uint32_t k1 = k0 + 1;

				x = std::clamp(x, 0.5, fluid.N + 0.5);
				i0 = static_cast<std::uint32_t>(x);
				i1 = i0 + 1;
				y = std::clamp(y, 0.5, fluid.N + 0.5);
				j0 = static_cast<std::uint32_t>(y);
				j1 = j0 + 1;
				z = std::clamp(z, 0.5, fluid.N + 0.5);
				k0 = static_cast<std::uint32_t>(z);
				k1 = k0 + 1;

				double s1 = x	- i0;
				double s0 = 1.0 - s1;
				double t1 = y	- j0;
				double t0 = 1.0 - t1;
				double u1 = z	- k0;
				double u0 = 1.0 - u1;


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
    // Reverse advection to calculate errors made, than correct the first advection to reduce the errors
	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t j = 1; j <= fluid.N; ++j)
		{
			for (std::uint32_t i = 1; i <= fluid.N; ++i)
			{
				double x = i-dt0*X[fluid.IX(i,j,k)];
				double y = j-dt0*Y[fluid.IX(i,j,k)];
				double z = k-dt0*Z[fluid.IX(i,j,k)];
				x = std::clamp(x, 0.5, fluid.N + 0.5);
				y = std::clamp(y, 0.5, fluid.N + 0.5);
				z = std::clamp(z, 0.5, fluid.N + 0.5);

				x = std::clamp(x, 0.5, fluid.N + 0.5);
				std::uint32_t i0 = static_cast<std::uint32_t>(x);
				std::uint32_t i1 = i0 + 1;
				y = std::clamp(y, 0.5, fluid.N + 0.5);
				std::uint32_t j0 = static_cast<std::uint32_t>(y);
				std::uint32_t j1 = j0 + 1;
				z = std::clamp(z, 0.5, fluid.N + 0.5);
				std::uint32_t k0 = static_cast<std::uint32_t>(z);
				std::uint32_t k1 = k0 + 1;

				x = std::clamp(x, 0.5, fluid.N + 0.5);
				i0 = static_cast<std::uint32_t>(x);
				i1 = i0 + 1;
				y = std::clamp(y, 0.5, fluid.N + 0.5);
				j0 = static_cast<std::uint32_t>(y);
				j1 = j0 + 1;
				z = std::clamp(z, 0.5, fluid.N + 0.5);
				k0 = static_cast<std::uint32_t>(z);
				k1 = k0 + 1;

                double top = std::max({D[fluid.IX(i0,j0,k0)],D[fluid.IX(i0,j0,k1)],D[fluid.IX(i0,j1,k0)],D[fluid.IX(i0,j1,k1)],D[fluid.IX(i1,j0,k0)],D[fluid.IX(i1,j0,k1)],D[fluid.IX(i1,j1,k0)],D[fluid.IX(i1,j1,k1)]});
                double bot = std::min({D[fluid.IX(i0,j0,k0)],D[fluid.IX(i0,j0,k1)],D[fluid.IX(i0,j1,k0)],D[fluid.IX(i0,j1,k1)],D[fluid.IX(i1,j0,k0)],D[fluid.IX(i1,j0,k1)],D[fluid.IX(i1,j1,k0)],D[fluid.IX(i1,j1,k1)]});

				x = i+dt0*X[fluid.IX(i,j,k)];
				y = j+dt0*Y[fluid.IX(i,j,k)];
				z = k+dt0*Z[fluid.IX(i,j,k)];
				x = std::clamp(x, 0.5, fluid.N + 0.5);
				y = std::clamp(y, 0.5, fluid.N + 0.5);
				z = std::clamp(z, 0.5, fluid.N + 0.5);

				x = std::clamp(x, 0.5, fluid.N + 0.5);
				i0 = static_cast<std::uint32_t>(x);
				i1 = i0 + 1;
				y = std::clamp(y, 0.5, fluid.N + 0.5);
				j0 = static_cast<std::uint32_t>(y);
				j1 = j0 + 1;
				z = std::clamp(z, 0.5, fluid.N + 0.5);
				k0 = static_cast<std::uint32_t>(z);
				k1 = k0 + 1;

				x = std::clamp(x, 0.5, fluid.N + 0.5);
				i0 = static_cast<std::uint32_t>(x);
				i1 = i0 + 1;
				y = std::clamp(y, 0.5, fluid.N + 0.5);
				j0 = static_cast<std::uint32_t>(y);
				j1 = j0 + 1;
				z = std::clamp(z, 0.5, fluid.N + 0.5);
				k0 = static_cast<std::uint32_t>(z);
				k1 = k0 + 1;

				double s1 = x	- i0;
				double s0 = 1.0 - s1;
				double t1 = y	- j0;
				double t0 = 1.0 - t1;
				double u1 = z	- k0;
				double u0 = 1.0 - u1;

                double back = 
					s0*(t0*(u0*D[fluid.IX(i0,j0,k0)]
							+u1*D[fluid.IX(i0,j0,k1)])
						+(t1*(u0*D[fluid.IX(i0,j1,k0)]
							+u1*D[fluid.IX(i0,j1,k1)])))
					+s1*(t0*(u0*D[fluid.IX(i1,j0,k0)]
							+u1*D[fluid.IX(i1,j0,k1)])
						+(t1*(u0*D[fluid.IX(i1,j1,k0)]
							+u1*D[fluid.IX(i1,j1,k1)])));
                D[fluid.IX(i,j,k)] = std::clamp(D[fluid.IX(i,j,k)] + 0.5 * (Dprev[fluid.IX(i,j,k)] - back), bot, top);
			}
		}
	}
	setBnd(fluid, D, b);
}

void Fluids::GaussSeidelRelaxationLinSolve(const Fluid3D& fluid, std::vector<double>& X, std::vector<double>& Xprev, float a, float c, std::uint8_t b) const
{
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

void Fluids::applyPreconditioner(const Fluid3D& fluid, const Eigen::VectorXd& r, const Laplacian& A, Eigen::VectorXd& z) const
{
    // Solve Lq = r
    Eigen::VectorXd q = Eigen::VectorXd::Zero(z.size());
    for (std::uint32_t k = 0; k < fluid.N; ++k)
    {
		for (std::uint32_t j = 0; j < fluid.N; ++j)
        {
			for (std::uint32_t i = 0; i < fluid.N; ++i)
            {
                std::uint32_t ind   = i+j*fluid.N+k*fluid.N2;
                if (ind < z.size())
                {
                    std::uint32_t indmi = (i-1)+j*fluid.N+k*fluid.N2;
                    std::uint32_t indmj = i+(j-1)*fluid.N+k*fluid.N2;
                    std::uint32_t indmk = i+j*fluid.N+(k-1)*fluid.N2;

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
    for (std::int32_t k = fluid.N-1; k >= 0; --k)
    {
		for (std::int32_t j = fluid.N-1; j >= 0; --j)
        {
			for (std::int32_t i = fluid.N-1; i >= 0; --i)
            {
                std::uint32_t ind   = i+j*fluid.N+k*fluid.N2;
                if (ind < z.size())
                {
                    std::uint32_t indpi = (i+1)+j*fluid.N+k*fluid.N2;
                    std::uint32_t indpj = i+(j+1)*fluid.N+k*fluid.N2;
                    std::uint32_t indpk = i+j*fluid.N+(k+1)*fluid.N2;

                    double a = indpi < fluid.N3 ? z.coeff(indpi) : 0;
                    double b = indpj < fluid.N3 ? z.coeff(indpj) : 0;
                    double c = indpk < fluid.N3 ? z.coeff(indpk) : 0;

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
    std::uint32_t n = A.diag.size();
    Eigen::VectorXd x(n);
    Eigen::VectorXd b(n);

    // Filling matrices
    std::uint32_t bIt = 0;
    for (std::uint32_t k = 1; k <= fluid.N; ++k)
    {
        for (std::uint32_t j = 1; j <= fluid.N; ++j)
        {
            for (std::uint32_t i = 1; i <= fluid.N; ++i)
            {
                if (bIt >= n)
                {
                    break;
                }
                b.coeffRef(bIt) = Xprev[fluid.IX(i,j,k)];
                bIt++;
            }
        }
    }

    // Solving Ap = b
    Eigen::VectorXd r = b;
    if (r.isZero(0))
    {
        X = Xprev;
        return;
    }
    Eigen::VectorXd p = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd z = p;
    applyPreconditioner(fluid, r, A, z);
    Eigen::VectorXd s = z;
    double sig = z.dot(r);

    for (std::uint32_t i = 0; i < b.size(); ++i)
    {
        z = A.A * s;
        double alpha = sig / s.dot(z);
        p = p + alpha * s;
        r = r - alpha * z;
        if (r.lpNorm<Eigen::Infinity>() < 10e-5)
        {
            break;
        }
        applyPreconditioner(fluid, r, A, z);
        double signew = z.dot(r);
        double beta = signew / sig;
        s = z + beta * s;
        sig = signew;
    }

    // Write the results
    for (std::uint32_t k = 1; k <= fluid.N; ++k)
    {
        for (std::uint32_t j = 1; j <= fluid.N; ++j)
        {
            for (std::uint32_t i = 1; i <= fluid.N; ++i)
            {
                std::uint32_t ind = (i-1)+(j-1)*fluid.N+(k-1)*fluid.N*fluid.N;
                if (ind >= n)
                {
                    break;
                }
                X[fluid.IX(i,j,k)] = p.coeffRef(ind);
            }
        }
    }

    setBnd(fluid, X, bs);
}

void Fluids::project(const Fluid3D& fluid, std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Z, std::vector<double>& p, std::vector<double>& div)
{
    double h = 1.0/fluid.N;
	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t j = 1; j <= fluid.N; ++j)
		{
			for (std::uint32_t i = 1; i <= fluid.N; ++i)
			{
                div[fluid.IX(i,j,k)] = -(0.33)*
                    ((X[fluid.IX(i+1,j,k)]-X[fluid.IX(i-1,j,k)])*h+
                     (Y[fluid.IX(i,j+1,k)]-Y[fluid.IX(i,j-1,k)])*h+
                     (Z[fluid.IX(i,j,k+1)]-Z[fluid.IX(i,j,k-1)])*h);
				p[fluid.IX(i,j,k)] = 0;
			}
		}
	}
	setBnd(fluid, div, 0);
	setBnd(fluid, p, 0);

    ConjugateGradientMethodLinSolve(fluid, p, div, 0, fluid.laplacianProject);
    //GaussSeidelRelaxationLinSolve(fluid, p, div, 1, 6, 0);

	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t j = 1; j <= fluid.N; ++j)
		{
			for (std::uint32_t i = 1; i <= fluid.N; ++i)
			{
				X[fluid.IX(i,j,k)] -= 0.5*fluid.N*(p[fluid.IX(i+1,j,k)]-p[fluid.IX(i-1,j,k)]);
				Y[fluid.IX(i,j,k)] -= 0.5*fluid.N*(p[fluid.IX(i,j+1,k)]-p[fluid.IX(i,j-1,k)]);
				Z[fluid.IX(i,j,k)] -= 0.5*fluid.N*(p[fluid.IX(i,j,k+1)]-p[fluid.IX(i,j,k-1)]);
			}
		}
	}

	setBnd(fluid, X, 1);
	setBnd(fluid, Y, 2);
	setBnd(fluid, Z, 3);
}

void Fluids::setBnd(const Fluid3D& fluid, std::vector<double>& X, const std::uint8_t b) const
{
	// Faces
	for (std::uint32_t j = 1; j <= fluid.N; ++j)
	{
		for (std::uint32_t i = 1; i <= fluid.N; ++i)
		{
			X[fluid.IX(i,j,0)]		    = b == 3 ? -X[fluid.IX(i,j,1)]		    : X[fluid.IX(i,j,1)];
			X[fluid.IX(i,j,fluid.N)]    = b == 3 ? -X[fluid.IX(i,j,fluid.N-1)]	: X[fluid.IX(i,j,fluid.N-1)];
		}
	}
	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t i = 1; i <= fluid.N; ++i)
		{
			X[fluid.IX(i,0,k)]		    = b == 2 ? -X[fluid.IX(i,1,k)]		    : X[fluid.IX(i,1,k)];
			X[fluid.IX(i,fluid.N,k)]    = b == 2 ? -X[fluid.IX(i,fluid.N-1,k)]	: X[fluid.IX(i,fluid.N-1,k)];
		}
	}
	for (std::uint32_t k = 1; k <= fluid.N; ++k)
	{
		for (std::uint32_t j = 1; j <= fluid.N; ++j)
		{
			X[fluid.IX(0,j,k)]		    = b == 1 ? -X[fluid.IX(1,j,k)]		    : X[fluid.IX(1,j,k)];
			X[fluid.IX(fluid.N,j,k)]    = b == 1 ? -X[fluid.IX(fluid.N-1,j,k)]	: X[fluid.IX(fluid.N-1,j,k)];
		}
	}

	// Edges
	for (std::uint32_t i = 1; i <= fluid.N; ++i)
	{
		X[fluid.IX(i,0,0)]              = 0.5*(X[fluid.IX(i,1,0)]              +X[fluid.IX(i,0,1)]);
		X[fluid.IX(i,fluid.N,0)]        = 0.5*(X[fluid.IX(i,fluid.N-1,0)]      +X[fluid.IX(i,fluid.N,1)]);
		X[fluid.IX(i,0,fluid.N)]        = 0.5*(X[fluid.IX(i,0,fluid.N-1)]      +X[fluid.IX(i,1,fluid.N)]);
		X[fluid.IX(i,fluid.N,fluid.N)]  = 0.5*(X[fluid.IX(i,fluid.N-1,fluid.N)]+X[fluid.IX(i,fluid.N,fluid.N-1)]);

		X[fluid.IX(0,i,0)]              = 0.5*(X[fluid.IX(1,i,0)]              +X[fluid.IX(0,i,1)]);
		X[fluid.IX(fluid.N,i,0)]        = 0.5*(X[fluid.IX(fluid.N-1,i,0)]      +X[fluid.IX(fluid.N,i,1)]);
		X[fluid.IX(0,i,fluid.N)]        = 0.5*(X[fluid.IX(0,i,fluid.N-1)]      +X[fluid.IX(1,i,fluid.N)]);
		X[fluid.IX(fluid.N,i,fluid.N)]  = 0.5*(X[fluid.IX(fluid.N-1,i,fluid.N)]+X[fluid.IX(fluid.N,i,0)]);

		X[fluid.IX(0,0,i)]              = 0.5*(X[fluid.IX(0,1,i)]              +X[fluid.IX(1,0,i)]);
		X[fluid.IX(0,fluid.N,i)]        = 0.5*(X[fluid.IX(0,fluid.N-1,i)]      +X[fluid.IX(1,fluid.N,i)]);
		X[fluid.IX(fluid.N,0,i)]        = 0.5*(X[fluid.IX(fluid.N-1,0,i)]      +X[fluid.IX(fluid.N,1,i)]);
		X[fluid.IX(fluid.N,fluid.N,i)]  = 0.5*(X[fluid.IX(fluid.N,fluid.N-1,i)]+X[fluid.IX(fluid.N-1,fluid.N,i)]);
	}

	// Corners
	X[fluid.IX(0,0,0)]				    = 0.33*(X[fluid.IX(1,0,0)]             +X[fluid.IX(0,1,0)]                 +X[fluid.IX(0,0,1)]);
	X[fluid.IX(0,fluid.N,0)]		    = 0.33*(X[fluid.IX(1,fluid.N,0)]       +X[fluid.IX(0,fluid.N-1,0)]         +X[fluid.IX(0,fluid.N,1)]);
	X[fluid.IX(0,0,fluid.N)]		    = 0.33*(X[fluid.IX(1,0,fluid.N)]       +X[fluid.IX(0,1,fluid.N)]           +X[fluid.IX(0,0,fluid.N+1)]);
	X[fluid.IX(0,fluid.N,fluid.N)]      = 0.33*(X[fluid.IX(1,fluid.N,fluid.N)] +X[fluid.IX(0,fluid.N-1,fluid.N)]   +X[fluid.IX(0,fluid.N,fluid.N-1)]);

	X[fluid.IX(fluid.N,0,0)]				= 0.33*(X[fluid.IX(fluid.N-1,0,0)]             +X[fluid.IX(fluid.N,1,0)]               +X[fluid.IX(fluid.N,0,1)]);
	X[fluid.IX(fluid.N,fluid.N,0)]		    = 0.33*(X[fluid.IX(fluid.N-1,fluid.N,0)]       +X[fluid.IX(fluid.N,fluid.N-1,0)]       +X[fluid.IX(fluid.N,fluid.N,1)]);
	X[fluid.IX(fluid.N,0,fluid.N)]		    = 0.33*(X[fluid.IX(fluid.N-1,0,fluid.N)]       +X[fluid.IX(fluid.N,1,fluid.N)]         +X[fluid.IX(fluid.N,0,fluid.N-1)]);
	X[fluid.IX(fluid.N,fluid.N,fluid.N)]    = 0.33*(X[fluid.IX(fluid.N-1,fluid.N,fluid.N)] +X[fluid.IX(fluid.N,fluid.N-1,fluid.N)] +X[fluid.IX(fluid.N,fluid.N,fluid.N-1)]);
}

void Fluids::updateRender(Fluid3D& fluid)
{
	std::vector<std::uint8_t> texture((fluid.N+2)*(fluid.N+2)*(fluid.N+2), 0);

	for (std::uint32_t i = 0; i < (fluid.N+2)*(fluid.N+2)*(fluid.N); ++i)
	{
		texture[i] = static_cast<std::uint8_t>(std::clamp(fluid.substanceField[i], 0.0, 255.0));
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
        debugViscosity = fluid.viscosity;
        debugDiffusion = fluid.diffusion;
        debugDt = fluid.dt;
        debugN = fluid.N;

		auto& material = gCoordinator.GetComponent<Material>(entity);
        debugLightIntensity[0] = material.lightIntensity.x;
        debugLightIntensity[1] = material.lightIntensity.y;
        debugLightIntensity[2] = material.lightIntensity.z;
        debugAbsorption = material.absorption;
    }
}

void Fluids::fluidDebugTool()
{
    ImGui::Begin("Fluid");
    ImGui::SliderInt("N", &debugN, 4, 128);
    ImGui::SliderFloat("viscosity", &debugViscosity, 0.0f, 16.0f);
    ImGui::SliderFloat("diffusion", &debugDiffusion, 0.0f, 16.0f);
    ImGui::SliderFloat("dt", &debugDt, 0.000001f, 0.0002f, "%.8f");
    ImGui::SliderFloat("absorption", &debugAbsorption, 1.0f, 100.0f);
    ImGui::SliderFloat3("lightIntensity", debugLightIntensity, 0.0f, 1.0f);
    if (ImGui::Button("Apply"))
    {
        reset();
    }
    if (ImGui::Button("Reset"))
    {
        reset(true);
    }
    ImGui::End();
}
#endif
