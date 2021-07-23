#pragma once

#include "types.h"
#include "StaggeredGrid.h"
#include "ConjugateGradient.h"

#include <chrono>
#include <unordered_map>
#include <algorithm>
#include <execution>

struct Cell
{
    std::uint16_t i;
    std::uint16_t j;
    std::uint64_t label;
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

    inline void advect(Field<double,std::uint16_t>& F, const Field<double,std::uint16_t>& Fprev, const std::uint8_t b) const;
    void project(Field<double,std::uint16_t>& X, Field<double,std::uint16_t>& Y);

    void setBnd(Field<double,std::uint16_t>& F, const std::uint8_t b) const;

    void updateTexture();

    void reinitLevelSet(const std::uint64_t nbIte);
    double gradLength(const Field<double,std::uint16_t>& F, const std::uint16_t i, const std::uint16_t j) const;

    void writeVolumeFile(const std::uint64_t iteration);

    void extrapolate(Field<double,std::uint16_t>& F);

    inline std::uint64_t hash(const std::uint16_t i, const std::uint16_t j) const;

    void initCG();

    inline double interp(const Field<double,std::uint16_t>& F, double x, double y) const;

    void setActiveCells();

    std::vector<glm::vec2> particles {};

    constexpr static const std::uint16_t _N = 65;
    constexpr static const double _viscosity = 1.15;
    constexpr static const double _diffusion = 0.000;
    constexpr static const double _dt = 0.0005;
    constexpr static const Solver _solverType = CG;
    constexpr static const Advection _advectionType = SEMI_LAGRANGIAN;

    Laplacian _laplacianProject {};
    Laplacian _laplacianViscosityX {};
    Laplacian _laplacianViscosityY {};
    Laplacian _laplacianDiffuse {};

    std::vector<std::uint8_t> _texture = std::vector<std::uint8_t>(_N*_N*3);

    StaggeredGrid<double, std::uint16_t> _grid {_N};

    std::unordered_map<std::uint64_t, Cell> _activeCells;
};

template<typename T>
void write(std::ofstream &f, const T data)
{
    f.write(reinterpret_cast<const char *>(&data), sizeof(data));
}
