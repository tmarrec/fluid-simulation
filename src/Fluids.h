#pragma once

#include <chrono>
#include <unordered_map>
#include <algorithm>
#include <execution>
#include <vector>
#include <memory>

#include "./types.h"
#include "./config.h"
#include "./StaggeredGrid.h"
#include "./Advect.h"
#include "./Project.h"

class Fluids
{
 public:
    explicit Fluids();
    void update(const std::uint64_t iteration);
    const std::vector<std::uint8_t>& texture() const;
    const std::vector<double>& X() const;
    const std::vector<double>& Y() const;
    const Field<double, std::uint16_t>& surface() const;
    bool isCellActive(
            const std::uint16_t i,
            const std::uint16_t j,
            const std::uint16_t k
        ) const;

 private:
    void step();
    void addForces();
    void redistancing(
            const std::uint64_t nbIte,
            Field<double, std::uint16_t>& field,
            Field<double, std::uint16_t>& fieldTemp
        ) const;
    void extrapolate(
            Field<double, std::uint16_t>& F,
            Field<double, std::uint16_t>& Ftemp,
            std::uint16_t nbIte = 0
        ) const;
    void updateTexture2D();
    void updateTexture3D();

    std::uint64_t _iteration = 0;
    std::vector<std::uint8_t> _texture;
    StaggeredGrid<double, std::uint16_t> _grid {Config::N};

    std::unique_ptr<Advect> _advection;
    std::unique_ptr<Project> _projection;
};

