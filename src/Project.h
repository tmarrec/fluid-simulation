#pragma once

#include "types.h"
#include "StaggeredGrid.h"
#include "ConjugateGradient.h"

class Project
{
public:
    explicit Project(StaggeredGrid<double, std::uint16_t>& grid) : _grid(grid) {};
    virtual void project() = 0;
protected:
    inline double pressureAt(const std::uint16_t i, const std::uint16_t j, const std::uint16_t k, const Eigen::VectorXd x, const std::uint64_t l) const;

    Eigen::SparseMatrix<double> _A;
    StaggeredGrid<double, std::uint16_t>& _grid;
};

class Project2D : public Project
{
using Project::Project;
public:
    virtual void project();
private:
    void preparePressureSolving(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b);
    inline double div(const std::uint16_t i, const std::uint16_t j, const std::uint16_t k) const;
};


class Project3D : public Project
{
using Project::Project;
public:
    virtual void project();
private:
    void preparePressureSolving(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b);
    inline double div(const std::uint16_t i, const std::uint16_t j, const std::uint16_t k) const;
};
