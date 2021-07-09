#pragma once

#include "../ecs/Coordinator.h"
#include "../Components.h"
#include "../BasicEntities.h"

#include <Eigen/Sparse>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/IterativeLinearSolvers/BasicPreconditioners.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <array>
#include <cstdint>
#include <vector>
#include <unistd.h>
#include <ctime>
#include <numeric>
#include <iomanip>

extern Coordinator gCoordinator;

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

template<typename T, typename U>
class Field
{
public:
    Field(U Xsize, U Ysize) : _Xsize(Xsize), _Ysize(Ysize)
    {
        _grid.resize(Xsize*Ysize);
        _maxIt = Xsize*Ysize;
    }
            T& operator()(U idx)            { return _grid[idx]; }
    const   T& operator()(U idx)    const   { return _grid[idx]; }
            T& operator()(U i, U j)         { return _grid[_idx(i,j)]; }
    const   T& operator()(U i, U j) const   { return _grid[_idx(i,j)]; }
    const   U& x()                  const   { return _Xsize; }
    const   U& y()                  const   { return _Ysize; }
    const   std::uint64_t& maxIt()  const   { return _maxIt; }
    const   std::vector<T>& data()  const   { return _grid; }

    friend std::ostream& operator<<(std::ostream& os, const Field& obj)
    {
        os << std::fixed << std::setprecision(2);
        for (U j = 0; j < obj._Ysize; ++j)
        {
            for (U i = 0; i < obj._Xsize; ++i)
            {
                if (obj(i,j) >= 0)
                {
                    os << " ";
                }
                os << obj(i,j) << " ";
            }
            os << std::endl;
        }
        return os;
    }

private:
    inline std::uint64_t _idx(const U i, const U j) const
    { 
        return i + j * _Xsize;
    };

    std::vector<T> _grid;
    U _Xsize;
    U _Ysize;
    std::uint64_t _maxIt;
};

class Fluids : public System
{
public:
    Fluids();
    void init(std::shared_ptr<Renderer> renderer);
    void update([[maybe_unused]] std::uint64_t iteration);

private:
    void Vstep();
    void Sstep();

    void addSource(std::vector<double>& X, const std::vector<double>& S) const;
    void diffuse(Field<double,std::uint16_t>& F, const Field<double,std::uint16_t>& Fprev, const std::uint8_t b, const Laplacian& A);

    void advect(Field<double,std::uint16_t>& F, const Field<double,std::uint16_t>& Fprev, const Field<double,std::uint16_t>& X, const Field<double,std::uint16_t>& Y, const std::uint8_t b) const;
    void project(Field<double,std::uint16_t>& X, Field<double,std::uint16_t>& Y, Field<double,std::uint16_t>& p, Field<double,std::uint16_t>& div);

    //void GaussSeidelRelaxationLinSolve(const Fluid3D& fluid, std::vector<double>& X, const std::vector<double>& Xprev, const double a, const double c, std::uint8_t b) const;
    void ConjugateGradientMethodLinSolve(Field<double,std::uint16_t>& F, const Field<double,std::uint16_t>& Fprev, const std::uint8_t bs, const Laplacian& A);

    void applyPreconditioner(const Eigen::VectorXd& r, const Laplacian& A, Eigen::VectorXd& z, const Solver solver) const;

    void setBnd(Field<double,std::uint16_t>& F, const std::uint8_t b) const;

    void updateRender(Fluid3D& fluid);

    /*
    double gradLength(const Fluid3D& fluid, const std::vector<double>& X, const std::uint64_t i, const std::uint64_t j) const;
    void reinitLevelSet(Fluid3D& fluid, const std::uint64_t nbIte) const;
    */

    void writeVolumeFile(std::uint64_t iteration);


    void setAMatrices(Laplacian& laplacian) const;

    void setPrecon(Laplacian& A) const;

    void initCG();


    std::shared_ptr<Renderer> _renderer = nullptr;

    double VstepTime = 0;
    double SstepTime = 0;
    double VstepProjectTime = 0;
    double VstepAdvectTime = 0;
    double VstepDiffuseTime = 0;
    double SstepDiffuseTime = 0;
    double SstepAdvectTime = 0;
    std::vector<glm::vec2> particles;

    constexpr static const std::uint16_t _N = 513;
    constexpr static const double _viscosity = 1.15;
    constexpr static const double _diffusion = 0.000;
    constexpr static const double _dt = 0.0005;
    constexpr static const Solver _solverType = CG;
    constexpr static const Advection _advectionType = SEMI_LAGRANGIAN;

    Field<double, std::uint16_t> _substance {_N, _N};
    Field<double, std::uint16_t> _prevSubstance {_N, _N};
    Field<double, std::uint16_t> _fieldX {_N+1, _N};
    Field<double, std::uint16_t> _fieldY {_N, _N+1};
    Field<double, std::uint16_t> _prevFieldX {_N+1, _N};
    Field<double, std::uint16_t> _prevFieldY {_N, _N+1};

    Laplacian _laplacianProject {};
    Laplacian _laplacianViscosityX {};
    Laplacian _laplacianViscosityY {};
    Laplacian _laplacianDiffuse {};
};

template<typename T>
void write(std::ofstream &f, T data)
{
    f.write(reinterpret_cast<const char *>(&data), sizeof(data));
}
