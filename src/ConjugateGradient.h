#pragma once

#include "types.h"

#include "StaggeredGrid.h"

#include <Eigen/Sparse>

struct Laplacian
{
    Eigen::SparseMatrix<double> A; 
    //Eigen::VectorXd precon;
    Eigen::VectorXd indmi;
    Eigen::VectorXd indmj;
    Eigen::VectorXd indmk;
    Eigen::VectorXd indpi;
    Eigen::VectorXd indpj;
    Eigen::VectorXd indpk;
    Eigen::VectorXd i;
    Eigen::VectorXd j;
    Eigen::VectorXd k;
    std::uint8_t minus = 0;

    /*
    Eigen::VectorXd Adiag;
    Eigen::VectorXd Ax;
    Eigen::VectorXd Ay;
    Eigen::VectorXd Az;
    std::uint16_t n;
    */
};

void ConjugateGradient(const Laplacian& A, Eigen::VectorXd& x, const Eigen::VectorXd& b, Solver solverType, StaggeredGrid<double, std::uint16_t>& grid);
void applyPreconditioner(const Eigen::VectorXd& r, Eigen::VectorXd& z, const Solver solver, StaggeredGrid<double, std::uint16_t>& grid);
void buildPrecondtioner(StaggeredGrid<double, std::uint16_t>& grid);
