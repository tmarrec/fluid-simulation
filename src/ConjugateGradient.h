#pragma once

#include "types.h"

#include <Eigen/Sparse>

struct Laplacian
{
    Eigen::SparseMatrix<double> A; 
    Eigen::VectorXd diag;
    Eigen::VectorXd plusi;
    Eigen::VectorXd plusj;
    Eigen::VectorXd plusk;
    Eigen::VectorXd precon;
    std::uint8_t minus = 0;
};

void ConjugateGradient(const Laplacian& A, Eigen::VectorXd& x, const Eigen::VectorXd& b, Solver solverType);
void applyPreconditioner(const Eigen::VectorXd& r, const Laplacian& A, Eigen::VectorXd& z, const Solver solver);
void setPrecon(Laplacian& A);
void setAMatrices(Laplacian& laplacian);
void initLaplacians(const std::uint16_t N, const double dt, const double viscosity, const double diffusion, Laplacian& laplacianViscosityX, Laplacian& laplacianViscosityY, Laplacian& laplacianProject, Laplacian& laplacianDiffuse);
