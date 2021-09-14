#pragma once

#include "types.h"
#include "config.h"

#include "StaggeredGrid.h"

#include <Eigen/Sparse>

void ConjugateGradient(const Eigen::SparseMatrix<double>& A, Eigen::VectorXd& x, const Eigen::VectorXd& b, StaggeredGrid<double, std::uint16_t>& grid);
void applyPreconditioner(const Eigen::VectorXd& r, Eigen::VectorXd& z, StaggeredGrid<double, std::uint16_t>& grid);
void buildPrecondtioner(StaggeredGrid<double, std::uint16_t>& grid);
