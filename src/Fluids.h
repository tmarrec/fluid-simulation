#pragma once

#include "types.h"

#include <Eigen/Sparse>
#include <iomanip>
#include <chrono>

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
        _labels.resize(Xsize*Ysize);
        _gridSet.resize(Xsize*Ysize);
        for (std::uint64_t it = 0; it < Xsize*Ysize; ++it)
        {
            _grid[it] = 0;
            _labels[it] = 0;
            _gridSet[it] = false;
        }
        _maxIt = Xsize*Ysize;
    }
            T& operator()(const U idx)                                  { return _grid[idx]; }
    const   T& operator()(const U idx)                          const   { return _grid[idx]; }
            T& operator()(const U i, const U j)                         { return _grid[idx(i,j)]; }
    const   T& operator()(const U i, const U j)                 const   { return _grid[idx(i,j)]; }
    const   U& x()                                              const   { return _Xsize; }
    const   U& y()                                              const   { return _Ysize; }
    const   std::uint64_t& maxIt()                              const   { return _maxIt; }
    const   std::vector<T>& data()                              const   { return _grid; }
            bool isSet(const U i, const U j)                    const   { return _gridSet[idx(i,j)]; }
            void set(const U i, const U j, const bool value)            { _gridSet[idx(i,j)] = value; }
            void set(const U idx, const bool value)                     { _gridSet[idx] = value; }
            void reset()
            {
                std::fill(_grid.begin(), _grid.end(), 0.0);
                for (std::uint64_t it = 0; it < _maxIt; ++it)
                {
                    set(it, false);
                }
            };
            void resetBool()
            {
                for (std::uint64_t it = 0; it < _maxIt; ++it)
                {
                    set(it, false);
                }
            };
            std::uint64_t& label(const U i, const U j)                  { return _labels[idx(i,j)]; }
    const   std::uint64_t& label(const U i, const U j)          const   { return _labels[idx(i,j)]; }
            void resetLabels()                                          { std::fill(_labels.begin(), _labels.end(), 0); }
            void setFromVec(const Eigen::VectorXd& v)
            {
                for (std::uint64_t it = 0; it < _maxIt; ++it)
                {
                    _grid[it] = v.coeff(it);
                }
            }
            Eigen::VectorXd vec()                               const
            { 
                Eigen::VectorXd v(_maxIt);
                for (std::uint64_t it = 0; it < _maxIt; ++it)
                {
                    v.coeffRef(it) = _grid[it];
                }
                return v;
            }

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

    inline std::uint64_t idx(const U i, const U j) const
    { 
        return i + j * _Xsize;
    };

private:
    std::vector<T> _grid;
    std::vector<std::uint64_t> _labels;
    std::vector<bool> _gridSet;
    U _Xsize;
    U _Ysize;
    std::uint64_t _maxIt;
};

class Fluids
{
public:
    Fluids();
    void update([[maybe_unused]] std::uint64_t iteration);
    const std::vector<std::uint8_t>& texture() const;
    const std::vector<double>& X() const;
    const std::vector<double>& Y() const;
    const std::uint16_t& N() const;

private:
    void vStep();
    void sStep();
    void levelSetStep();

    void diffuse(Field<double,std::uint16_t>& F, const Field<double,std::uint16_t>& Fprev, const std::uint8_t b, const Laplacian& A);

    void advect(Field<double,std::uint16_t>& F, const Field<double,std::uint16_t>& Fprev, const Field<double,std::uint16_t>& X, const Field<double,std::uint16_t>& Y, const std::uint8_t b) const;
    void project(Field<double,std::uint16_t>& X, Field<double,std::uint16_t>& Y);

    void ConjugateGradient(const Laplacian& A, Eigen::VectorXd& x, const Eigen::VectorXd& b);

    void applyPreconditioner(const Eigen::VectorXd& r, const Laplacian& A, Eigen::VectorXd& z, const Solver solver) const;

    void setBnd(Field<double,std::uint16_t>& F, const std::uint8_t b) const;

    void updateTexture();

    void reinitLevelSet(const std::uint64_t nbIte);
    double gradLength(const Field<double,std::uint16_t>& F, const std::uint16_t i, const std::uint16_t j) const;

    void writeVolumeFile(const std::uint64_t iteration);

    void extrapolate(Field<double,std::uint16_t>& F);


    void setAMatrices(Laplacian& laplacian) const;

    void setPrecon(Laplacian& A) const;

    void initCG();

    inline double getPressure(const std::uint16_t i, const std::uint16_t j, const std::uint16_t i2, const std::uint16_t j2);

    double VstepTime = 0;
    double SstepTime = 0;
    double VstepProjectTime = 0;
    double VstepAdvectTime = 0;
    double VstepDiffuseTime = 0;
    double SstepDiffuseTime = 0;
    double SstepAdvectTime = 0;
    std::vector<glm::vec2> particles {};

    constexpr static const std::uint16_t _N = 129;
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
    Field<double, std::uint16_t> _implicit {_N, _N};
    Field<double, std::uint16_t> _prevImplicit {_N, _N};
    Field<double, std::uint16_t> _p {_N, _N};
    Field<double, std::uint16_t> _div {_N, _N};

    Laplacian _laplacianProject {};
    Laplacian _laplacianViscosityX {};
    Laplacian _laplacianViscosityY {};
    Laplacian _laplacianDiffuse {};

    std::vector<std::uint8_t> _texture = std::vector<std::uint8_t>(_N*_N*3);
};

template<typename T>
void write(std::ofstream &f, const T data)
{
    f.write(reinterpret_cast<const char *>(&data), sizeof(data));
}
