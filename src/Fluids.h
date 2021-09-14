#pragma once

#include "types.h"
#include "config.h"
#include "StaggeredGrid.h"
#include "Advect.h"
#include "Project.h"

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
    const std::vector<double>& Z() const;
    const Field<double,std::uint16_t>& surface() const;
    bool isCellActive(std::uint16_t i, std::uint16_t j, std::uint16_t k) const;

private:
    void step();
    void addForces();
    void diffuse(Field<double,std::uint16_t>& F, const Field<double,std::uint16_t>& Fprev, const std::uint8_t b, const Eigen::SparseMatrix<double>& A);
    void redistancing(const std::uint64_t nbIte, Field<double, std::uint16_t>& field, Field<double, std::uint16_t>& fieldTemp);
    void writeVolumeFile(const std::uint64_t iteration);
    void extrapolate(Field<double,std::uint16_t>& F, Field<double,std::uint16_t>& Ftemp, std::uint16_t nbIte = 0);
    void updateTexture2D();
    void updateTexture3D();

    std::uint64_t _iteration = 0;
    std::vector<std::uint8_t> _texture;
    StaggeredGrid<double, std::uint16_t> _grid {Config::N};

    std::unique_ptr<Advect> _advection;
    std::unique_ptr<Project> _projection;
};

template<typename T>
void write(std::ofstream &f, const T data)
{
    f.write(reinterpret_cast<const char *>(&data), sizeof(data));
}
