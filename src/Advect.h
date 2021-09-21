#pragma once

#include "./types.h"
#include "./config.h"
#include "./StaggeredGrid.h"

class Advect
{
 public:
    virtual void advect(
            const StaggeredGrid<double, std::uint16_t>& grid,
            Field<double, std::uint16_t>& F,
            Field<double, std::uint16_t>& Fprev,
            const std::uint8_t b
        ) = 0;
};

class Advect2D : public Advect
{
 public:
    virtual void advect(
            const StaggeredGrid<double, std::uint16_t>& grid,
            Field<double, std::uint16_t>& F,
            Field<double, std::uint16_t>& Fprev,
            const std::uint8_t b
        ) override;
 private:
    inline double interp(
            const Field<double, std::uint16_t>& F,
            const double x,
            const double y
        ) const;
};


class Advect3D : public Advect
{
 public:
    virtual void advect(
            const StaggeredGrid<double, std::uint16_t>& grid,
            Field<double, std::uint16_t>& F,
            Field<double, std::uint16_t>& Fprev,
            const std::uint8_t b
        ) override;
 private:
    inline double interp(
            const Field<double, std::uint16_t>& F,
            const double x,
            const double y,
            const double z
        ) const;
};
