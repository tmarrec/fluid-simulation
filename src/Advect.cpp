#include "Advect.h"

void Advect3D::advect(StaggeredGrid<double, std::uint16_t>& grid, Field<double,std::uint16_t>& F, Field<double,std::uint16_t>& Fprev, const std::uint8_t b)
{
    Fprev = F;
	const double dt = Config::dt * Config::N;
    #pragma omp parallel for
    for (std::uint64_t n = 0; n < grid._surface.maxIt(); ++n)
    {
        const std::uint64_t xy = grid._surface.x()*grid._surface.y();
        const std::uint64_t m = n % xy;
        const std::uint16_t i = m % grid._surface.x();
        const std::uint16_t j = m / grid._surface.x();
        const std::uint16_t k = n / xy;
        if (!(F.label(i,j,k) & SOLID))
        {
            const double x = std::clamp((static_cast<double>(i)-dt*grid.getU(i,j,k,b)), 0.0, static_cast<double>(F.x()-0.0));
            const double y = std::clamp((static_cast<double>(j)-dt*grid.getV(i,j,k,b)), 0.0, static_cast<double>(F.y()-0.0));
            const double z = std::clamp((static_cast<double>(k)-dt*grid.getW(i,j,k,b)), 0.0, static_cast<double>(F.z()-0.0));

            F(i,j,k) = interp(Fprev, x, y, z);
            /*
            if (j == F.y()-1 && F(i,j,k) < 10 && b == 0)
            {
                std::cout << F << std::endl;
                interp(Fprev, x, y, z);
                std::cout << F(i,j,k) << std::endl;
                std::cout << i << " " << j << " " << k << std::endl;
                std::cout << x << " " << y << " " << z << std::endl;
                exit(0);
            }
            */
        }
    }
}

inline double Advect3D::interp(const Field<double,std::uint16_t>& F, double x, double y, double z) const
{
    const std::uint16_t i0 = static_cast<std::uint16_t>(x);
    const std::uint16_t i1 = std::clamp(i0 + 1, 1, int(F.x()-1));
    const std::uint16_t j0 = static_cast<std::uint16_t>(y);
    const std::uint16_t j1 = std::clamp(j0 + 1, 1, int(F.y()-1));
    const std::uint16_t k0 = static_cast<std::uint16_t>(z);
    const std::uint16_t k1 = std::clamp(k0 + 1, 1, int(F.z()-1));

    //std::cout << "> "<< j1 << std::endl;

    const double s1 = x - i0;
    const double s0 = 1.0 - s1;
    const double t1 = y - j0;
    const double t0 = 1.0 - t1;
    const double u1 = z - k0;
    const double u0 = 1.0 - u1;

    /*
    std::cout << k1 << std::endl;
    std::cout << "==" << std::endl;
    std::cout << F(i0,j0,k0) << std::endl;
    std::cout << F(i0,j0,k1) << std::endl;
    std::cout << F(i0,j1,k0) << std::endl;
    std::cout << F(i0,j1,k1) << std::endl;
    std::cout << F(i1,j0,k0) << std::endl;
    std::cout << F(i1,j0,k1) << std::endl;
    std::cout << F(i1,j1,k0) << std::endl;
    std::cout << F(i1,j1,k1) << std::endl;
    std::cout << "==" << std::endl;
    */

    return s0 * (     t0 * ( u0 * F(i0,j0,k0) + u1 * F(i0,j0,k1) )
                    + t1 * ( u0 * F(i0,j1,k0) + u1 * F(i0,j1,k1) )
                )
         + s1 * (     t0 * ( u0 * F(i1,j0,k0) + u1 * F(i1,j0,k1) )
                    + t1 * ( u0 * F(i1,j1,k0) + u1 * F(i1,j1,k1) )
                );
}

void Advect2D::advect(StaggeredGrid<double, std::uint16_t>& grid, Field<double,std::uint16_t>& F, Field<double,std::uint16_t>& Fprev, const std::uint8_t b)
{
    Fprev = F;
	const double dt = Config::dt * Config::N;
    #pragma omp parallel for
    for (std::uint64_t n = 0; n < grid._surface.maxIt(); ++n)
    {
        const std::uint64_t xy = grid._surface.x()*grid._surface.y();
        const std::uint64_t m = n % xy;
        const std::uint16_t i = m % grid._surface.x();
        const std::uint16_t j = m / grid._surface.x();
        if (!(F.label(i,j,0) & SOLID))
        {
            const double x = std::clamp((static_cast<double>(i)-dt*grid.getU(i,j,0,b)), 0.0, static_cast<double>(F.x()-0.0));
            const double y = std::clamp((static_cast<double>(j)-dt*grid.getV(i,j,0,b)), 0.0, static_cast<double>(F.y()-0.0));

            F(i,j,0) = interp(Fprev, x, y);
        }
    }
}

inline double Advect2D::interp(const Field<double,std::uint16_t>& F, double x, double y) const
{
    const std::uint16_t i0 = static_cast<std::uint16_t>(x);
    const std::uint16_t i1 = i0 + 1;
    const std::uint16_t j0 = static_cast<std::uint16_t>(y);
    const std::uint16_t j1 = j0 + 1;

    const double s1 = x - i0;
    const double s0 = 1.0 - s1;
    const double t1 = y - j0;
    const double t0 = 1.0 - t1;

    return s0 * (     t0 * F(i0,j0,0)
                    + t1 * F(i0,j1,0)
                )
         + s1 * (     t0 * F(i1,j0,0)
                    + t1 * F(i1,j1,0)
                );
}
