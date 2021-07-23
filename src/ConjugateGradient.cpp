#include "ConjugateGradient.h"

void ConjugateGradient(const Laplacian& A, Eigen::VectorXd& x, const Eigen::VectorXd& b, Solver solverType)
{
    const std::uint64_t diagSize = A.diag.size();

#ifdef DEBUG
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> lltOfA(A.A);
    if(lltOfA.info() == Eigen::NumericalIssue)
    {
        ERROR("DEBUG: Numerical Issue on the A matrix");
    }
#endif

    // Solving Ap = b
    Eigen::VectorXd r = b;
    if (r.isZero(0))
    {
        x = b;
        return;
    }
    x = Eigen::VectorXd::Zero(diagSize);
    
    Eigen::VectorXd z = x;
    applyPreconditioner(r, A, z, solverType);
    Eigen::VectorXd s = z;
    double sig = z.dot(r);

    for (std::uint64_t i = 0; i < static_cast<std::uint64_t>(b.size()); ++i)
    {
        z = A.A * s;
        const double alpha = sig / s.dot(z);
        x = x + alpha * s;
        r = r - alpha * z;
        if (r.lpNorm<Eigen::Infinity>() < 10e-5)
        {
            break;
        }
        applyPreconditioner(r, A, z, solverType);
        const double signew = z.dot(r);
        const double beta = signew / sig;
        s = z + beta * s;
        sig = signew;
    }
}

void applyPreconditioner(const Eigen::VectorXd& r, const Laplacian& A, Eigen::VectorXd& z, const Solver solver)
{
    if (solver == CG)
    {
        z = r;
        return;
    }
    ERROR("Preconditioner need rework");
    /*
    const std::uint64_t N2 = _N*_N; 
    const std::uint64_t N3 = _N*_N;

    // Solve Lq = r
    Eigen::VectorXd q = Eigen::VectorXd::Zero(z.size());
    for (std::int64_t n = 0; n < z.size(); ++n)
    {
        const std::uint64_t m = n % N2;
        const std::uint64_t i = m % _N;
        const std::uint64_t j = m / _N;
        const std::uint64_t k = n / N2;

        const std::uint64_t indmi = (i-1)+j*_N+k*N2;
        const std::uint64_t indmj = i+(j-1)*_N+k*N2;
        const std::uint64_t indmk = i+j*_N+(k-1)*N2;

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
        const std::uint64_t i = m % _N;
        const std::uint64_t j = m / _N;
        const std::uint64_t k = n / N2;

        const std::uint64_t indpi = (i+1)+j*_N+k*N2;
        const std::uint64_t indpj = i+(j+1)*_N+k*N2;
        const std::uint64_t indpk = i+j*_N+(k+1)*N2;

        const double a = indpi < N3 ? z.coeff(indpi) : 0;
        const double b = indpj < N3 ? z.coeff(indpj) : 0;
        const double c = indpk < N3 ? z.coeff(indpk) : 0;

        const double prec = A.precon.coeff(n);
        const double t = q.coeff(n) - A.plusi.coeff(n) * prec * a
                                    - A.plusj.coeff(n) * prec * b
                                    - A.plusk.coeff(n) * prec * c;
        z.coeffRef(n) = t * prec;
    }
    */
}


void setPrecon(Laplacian& A)
{
    /*
    const std::uint64_t N3 = A.A.rows();
    const std::uint64_t N2 = A.A.rows(); // TEMP
    A.precon = Eigen::VectorXd::Zero(N3);
    #pragma omp parallel for
    for (std::uint32_t n = 0; n < N3; ++n)
    {
        const std::uint32_t m = n % N2;
        const std::uint16_t i = m % _N;
        const std::uint16_t j = m / _N;
        const std::uint16_t k = n / N2;
        const std::uint32_t indmi = (i-1)+j*_N+k*N2;
        const std::uint32_t indmj = i+(j-1)*_N+k*N2;
        const std::uint32_t indmk = i+j*_N+(k-1)*N2;

        double a = 0;
        double b = 0;
        double c = 0;

        double i0 = 0;
        double i1 = 0;
        double i2 = 0;
        double i3 = 0;
        double j0 = 0;
        double j1 = 0;
        double j2 = 0;
        double j3 = 0;
        double k0 = 0;
        double k1 = 0;
        double k2 = 0;
        double k3 = 0;

        if (i > 0)
        {
            a = std::pow(A.plusi.coeff(indmi) * A.precon.coeff(indmi), 2);
            i0 = A.plusi.coeff(indmi);
            i1 = A.plusj.coeff(indmi);
            i2 = A.plusk.coeff(indmi);
            i3 = std::pow(A.precon.coeff(indmi), 2);
        }
        if (j > 0)
        {
            b = std::pow(A.plusj.coeff(indmj) * A.precon.coeff(indmj), 2);
            j0 = A.plusj.coeff(indmj);
            j1 = A.plusi.coeff(indmj);
            j2 = A.plusk.coeff(indmj);
            j3 = std::pow(A.precon.coeff(indmj), 2);
        }
        if (k > 0)
        {
            c = std::pow(A.plusk.coeff(indmk) * A.precon.coeff(indmk), 2);
            k0 = A.plusk.coeff(indmk);
            k1 = A.plusi.coeff(indmk);
            k2 = A.plusj.coeff(indmk);
            k3 = std::pow(A.precon.coeff(indmk), 2);
        }

        double e = A.diag.coeff(n) - a - b - c 
            - 0.97 * (
                    i0 * (i1 + i2) * i3
                +   j0 * (j1 + j2) * j3
                +   k0 * (k1 + k2) * k3
            );

        if (e < 0.25 * A.diag.coeff(n))
        {
            e = A.diag.coeff(n);
        }
        A.precon.coeffRef(n) = 1/std::sqrt(e);
    }
    */
}

void setAMatrices(Laplacian& laplacian)
{
    const std::uint64_t N3 = laplacian.A.rows();
    const std::uint64_t N2 = laplacian.A.rows(); // TEMP
    laplacian.diag = Eigen::VectorXd::Zero(N3);
    laplacian.plusi = Eigen::VectorXd::Zero(N3);
    laplacian.plusj = Eigen::VectorXd::Zero(N3);
    laplacian.plusk = Eigen::VectorXd::Zero(N3);

    /*
    #pragma omp parallel for
    for (std::uint32_t n = 0; n < N3; ++n)
    {
        const std::uint32_t m = n % N2;
        const std::uint16_t i = m % _N;
        const std::uint16_t j = m / _N;
        const std::uint16_t k = n / N2;
        laplacian.diag.coeffRef(n) = laplacian.A.coeff(n, n);
        if (n+1 < N3 && i+1 < _N)
        {
            laplacian.plusi[n] = laplacian.A.coeff(n, (i+1)+j*_N+k*N2);
        }
        if (n+_N < N3)
        {
            laplacian.plusj[n] = laplacian.A.coeff(n, i+(j+1)*_N+k*N2);
        }
        if (n+N2 < N3)
        {
            laplacian.plusk[n] = laplacian.A.coeff(n, i+j*_N+(k+1)*N2);
        }
    }
    */
}

