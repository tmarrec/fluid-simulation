#pragma once

#include "types.h"
#include "config.h"
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
    void addForces();
    void diffuse(Field<double,std::uint16_t>& F, const Field<double,std::uint16_t>& Fprev, const std::uint8_t b, const Laplacian& A);
    void advect(Field<double,std::uint16_t>& F, Field<double,std::uint16_t>& Fprev, const std::uint8_t b) const;
    void project();
    void updateTexture();
    void redistancing(const std::uint64_t nbIte, Field<double, std::uint16_t>& field, Field<double, std::uint16_t>& fieldTemp);
    void writeVolumeFile(const std::uint64_t iteration);
    void extrapolate(Field<double,std::uint16_t>& F, Field<double,std::uint16_t>& Ftemp, std::uint16_t nbIte = 0);
    void preparePressureSolving(Laplacian& A, Eigen::VectorXd& b);
    inline double interp(const Field<double,std::uint16_t>& F, double x, double y) const;
    inline double pressureAt(const std::uint16_t i, const std::uint16_t j, const Eigen::VectorXd x, const std::uint64_t l) const;
    inline double div(const std::uint16_t i, const std::uint16_t j) const;

    std::uint64_t _iteration = 0;
    std::vector<std::uint8_t> _texture = std::vector<std::uint8_t>(Config::N*Config::N*3);
    StaggeredGrid<double, std::uint16_t> _grid {Config::N};
};

template<typename T>
void write(std::ofstream &f, const T data)
{
    f.write(reinterpret_cast<const char *>(&data), sizeof(data));
}
