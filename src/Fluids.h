#pragma once

#include "types.h"
#include "StaggeredGrid.h"
#include "ConjugateGradient.h"

#include <chrono>
#include <unordered_map>
#include <algorithm>
#include <execution>

class Fluids
{
public:
    explicit Fluids();
    void update([[maybe_unused]] std::uint64_t iteration);
    const std::vector<std::uint8_t>& texture() const;
    const std::vector<double>& X() const;
    const std::vector<double>& Y() const;
    bool isCellActive(std::uint16_t i, std::uint16_t j) const;
    const std::uint16_t& N() const;

private:
    void step();
    void vStep();
    void sStep();
    void levelSetStep();

    void diffuse(Field<double,std::uint16_t>& F, const Field<double,std::uint16_t>& Fprev, const std::uint8_t b, const Laplacian& A);

    inline void advect(Field<double,std::uint16_t>& F, Field<double,std::uint16_t>& Fprev, const std::uint8_t b) const;
    void project();

    void setBnd(Field<double,std::uint16_t>& F, const std::uint8_t b) const;

    void updateTexture();

    void redistancing(const std::uint64_t nbIte, Field<double, std::uint16_t>& field);

    void writeVolumeFile(const std::uint64_t iteration);

    void extrapolate(Field<double,std::uint16_t>& F, std::uint16_t nbIte = 0);

    void pressureMatrix(Laplacian& A, Eigen::VectorXd& b) const;

    inline double interp(const Field<double,std::uint16_t>& F, double x, double y) const;

    inline double pressureAt(const std::uint16_t i, const std::uint16_t j, const Eigen::VectorXd x, const std::uint64_t l) const;

    std::vector<glm::vec2> particles {};

    constexpr static const std::uint16_t _N = 128;
    constexpr static const double _viscosity = 1.15;
    constexpr static const double _diffusion = 1.000;
    constexpr static const double _dt = 0.0005;
    constexpr static const Solver _solverType = CG;
    constexpr static const Advection _advectionType = SEMI_LAGRANGIAN;

    Laplacian _laplacianProject {};
    Laplacian _laplacianViscosityX {};
    Laplacian _laplacianViscosityY {};
    Laplacian _laplacianDiffuse {};

    std::vector<std::uint8_t> _texture = std::vector<std::uint8_t>(_N*_N*3);

    StaggeredGrid<double, std::uint16_t> _grid {_N};

    std::uint64_t _iteration = 0;

};

template<typename T>
void write(std::ofstream &f, const T data)
{
    f.write(reinterpret_cast<const char *>(&data), sizeof(data));
}
