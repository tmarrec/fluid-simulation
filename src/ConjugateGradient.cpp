#include "ConjugateGradient.h"

// Use the Conjugate Gradient method to solve the Ax = b system
void ConjugateGradient(
        const Eigen::SparseMatrix<double>& A,
        Eigen::VectorXd& x,
        const Eigen::VectorXd& b,
        StaggeredGrid<double, std::uint16_t>& grid
    )
{
    const std::uint64_t diagSize = x.size();

#ifdef DEBUG
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> lltOfA(A);
    if(lltOfA.info() == Eigen::NumericalIssue)
    {
        ERROR("DEBUG: Numerical Issue on the A matrix");
    }
#endif

    if (Config::solver == PCG)
    {
        buildPrecondtioner(grid);
    }

    // Solving Ap = b
    Eigen::VectorXd r = b;
    if (r.isZero(0))
    {
        x = b;
        return;
    }
    x = Eigen::VectorXd::Zero(diagSize);

    Eigen::VectorXd z = x;
    applyPreconditioner(r, z, grid);
    Eigen::VectorXd s = z;
    double sig = z.dot(r);

    for (std::uint64_t i = 0; i < static_cast<std::uint64_t>(b.size()); ++i)
    {
        z = A * s;
        const double alpha = sig / s.dot(z);
        x = x + alpha * s;
        r = r - alpha * z;
        if (r.lpNorm<Eigen::Infinity>() < 10e-5)
        {
            break;
        }
        applyPreconditioner(r, z, grid);
        const double signew = z.dot(r);
        const double beta = signew / sig;
        s = z + beta * s;
        sig = signew;
    }
}

// Create the preconditioner in the "_precon" grid of "grid"
void buildPrecondtioner(StaggeredGrid<double, std::uint16_t>& grid)
{
    for (std::int16_t k = 0; k < grid._surface.z(); ++k)
    {
        for (std::int16_t j = 0; j < grid._surface.y(); ++j)
        {
            for (std::int16_t i = 0; i < grid._surface.x(); ++i)
            {
                if (grid._surface.label(i, j, k) & LIQUID)
                {
                    double  a = 0.0,  b = 0.0,  c = 0.0;
                    double i0 = 0.0, i1 = 0.0, i2 = 0.0, i3 = 0.0;
                    double j0 = 0.0, j1 = 0.0, j2 = 0.0, j3 = 0.0;
                    double k0 = 0.0, k1 = 0.0, k2 = 0.0, k3 = 0.0;
                    if (i > 0.0)
                    {
                        a = std::pow(
                                grid._Ax(i-1, j, k) * grid._precon(i-1, j, k),
                                2
                            );
                        i0 = grid._Ax(i-1, j, k);
                        i1 = grid._Ay(i-1, j, k);
                        i2 = grid._Az(i-1, j, k);
                        i3 = std::pow(grid._precon(i-1, j, k), 2);
                    }
                    if (j > 0.0)
                    {
                        b = std::pow(
                                grid._Ay(i, j-1, k) * grid._precon(i, j-1, k),
                                2
                            );
                        j0 = grid._Ay(i, j-1, k);
                        j1 = grid._Ax(i, j-1, k);
                        j2 = grid._Az(i, j-1, k);
                        j3 = std::pow(grid._precon(i, j-1, k), 2);
                    }
                    if (k > 0.0)
                    {
                        c = std::pow(
                                grid._Az(i, j, k-1) * grid._precon(i, j, k-1),
                                2
                            );
                        k0 = grid._Az(i, j, k-1);
                        k1 = grid._Ax(i, j, k-1);
                        k2 = grid._Ay(i, j, k-1);
                        k3 = std::pow(grid._precon(i, j, k-1), 2);
                    }

                    double e = grid._Adiag(i, j, k) - a - b - c
                        - 0.97 * (
                                i0 * (i1 + i2) * i3
                            +   j0 * (j1 + j2) * j3
                            +   k0 * (k1 + k2) * k3
                        );

                    if (e < 0.25 * grid._Adiag(i, j, k))
                    {
                        e = grid._Adiag(i, j, k);
                    }
                    grid._precon(i, j, k) = 1.0/std::sqrt(e);
                }
            }
        }
    }
}

// Apply the previously computed preconditioner
// to the z vector to speed up CG convergence
void applyPreconditioner(
        const Eigen::VectorXd& r,
        Eigen::VectorXd& z,
        StaggeredGrid<double, std::uint16_t>& grid
    )
{
    if (Config::solver == CG)
    {
        z = r;
        return;
    }
    grid._q.reset();
    grid._z.reset();
    for (std::int16_t k = 0; k < grid._surface.z(); ++k)
    {
        for (std::int16_t j = 0; j < grid._surface.y(); ++j)
        {
            for (std::int16_t i = 0; i < grid._surface.x(); ++i)
            {
                if (grid._pressureID(i, j, k) > 0)
                {
                    const double a = i > 0
                        ?   (grid._Ax(i-1, j, k) *
                            grid._precon(i-1, j, k) *
                            grid._q(i-1, j, k))
                        : 0.0;
                    const double b = j > 0
                        ?   (grid._Ay(i, j-1, k) *
                            grid._precon(i, j-1, k) *
                            grid._q(i, j-1, k))
                        : 0.0;
                    const double c = k > 0
                        ?   (grid._Az(i, j, k-1) *
                            grid._precon(i, j, k-1) *
                            grid._q(i, j, k-1))
                        : 0.0;
                    const double t =
                        r(grid._pressureID(i, j, k)-1) - a - b - c;
                    grid._q(i, j, k) = t * grid._precon(i, j, k);
                }
            }
        }
    }
    for (std::int16_t k = grid._surface.z()-1; k >= 0; --k)
    {
        for (std::int16_t j = grid._surface.y()-1; j >= 0; --j)
        {
            for (std::int16_t i = grid._surface.x()-1; i >= 0; --i)
            {
                if (grid._surface.label(i, j, k) & LIQUID)
                {
                    const double a =
                        grid._Ax(i, j, k) *
                        grid._precon(i, j, k) *
                        grid._z(i+1, j, k);
                    const double b =
                        grid._Ay(i, j, k) *
                        grid._precon(i, j, k) *
                        grid._z(i, j+1, k);
                    const double c =
                        grid._Az(i, j, k) *
                        grid._precon(i, j, k) *
                        grid._z(i, j, k+1);

                    const double t = grid._q(i, j, k) - a - b - c;
                    grid._z(i, j, k) = t * grid._precon(i, j, k);
                    z(grid._pressureID(i, j, k)-1) = t * grid._precon(i, j, k);
                }
            }
        }
    }
}
