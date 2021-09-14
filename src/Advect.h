#pragma once

#include "types.h"
#include "config.h"
#include "StaggeredGrid.h"

class Advect
{
public:
    virtual void advect(StaggeredGrid<double, std::uint16_t>& grid, Field<double,std::uint16_t>& F, Field<double,std::uint16_t>& Fprev, const std::uint8_t b) = 0;
};

class Advect2D : public Advect
{
public:
    virtual void advect(StaggeredGrid<double, std::uint16_t>& grid, Field<double,std::uint16_t>& F, Field<double,std::uint16_t>& Fprev, const std::uint8_t b);
private:
    inline double interp(const Field<double,std::uint16_t>& F, double x, double y) const;
};


class Advect3D : public Advect
{
public:
    virtual void advect(StaggeredGrid<double, std::uint16_t>& grid, Field<double,std::uint16_t>& F, Field<double,std::uint16_t>& Fprev, const std::uint8_t b);
private:
    inline double interp(const Field<double,std::uint16_t>& F, double x, double y, double z) const;
};
