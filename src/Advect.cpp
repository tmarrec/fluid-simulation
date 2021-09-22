#include "Advect.h"

// 3D semi-lagrangian advection, going backward in time to get new values
void Advect3D::advect(
        const StaggeredGrid<double, std::uint16_t>& grid,
        Field<double, std::uint16_t>& F,
        Field<double, std::uint16_t>& Fprev,
        const std::uint8_t b
    )
{
    Fprev = F;
    const double dt = Config::dt * Config::N;
    #pragma omp parallel for
    for (std::uint64_t n = 0; n < grid._surface.maxIt(); ++n)
    {
        const std::uint64_t xy =
            static_cast<std::uint64_t>(grid._surface.x()*grid._surface.y());
        const std::uint64_t m = n % xy;
        const std::uint16_t i = m % grid._surface.x();
        const std::uint16_t j = m / grid._surface.x();
        const std::uint16_t k = n / xy;
        if (!(F.label(i, j, k) & SOLID))
        {
            const double x =
                std::clamp((static_cast<double>(i)-dt*grid.getU(i, j, k, b)),
                            0.0,
                            static_cast<double>(F.x()));
            const double y =
                std::clamp((static_cast<double>(j)-dt*grid.getV(i, j, k, b)),
                            0.0,
                            static_cast<double>(F.y()));
            const double z =
                std::clamp((static_cast<double>(k)-dt*grid.getW(i, j, k, b)),
                            0.0,
                            static_cast<double>(F.z()));
            F(i, j, k) = interp(Fprev, x, y, z);
        }
    }

    if (Config::advection == MACCORMACK)
    {
        // Reverse advection to calculate errors made,
        // than correct the first advection to reduce the errors
        #pragma omp parallel for
        for (std::uint64_t n = 0; n < grid._surface.maxIt(); ++n)
        {
            const std::uint64_t xy =
                static_cast<std::uint64_t>(grid._surface.x()*grid._surface.y());
            const std::uint64_t m = n % xy;
            const std::uint16_t i = m % grid._surface.x();
            const std::uint16_t j = m / grid._surface.x();
            const std::uint16_t k = n / xy;
            if (!(F.label(i, j, k) & SOLID))
            {
                double x = std::clamp(
                        (static_cast<double>(i)-dt*grid.getU(i, j, k, b)),
                        0.0,
                        static_cast<double>(F.x()));
                double y = std::clamp(
                        (static_cast<double>(j)-dt*grid.getV(i, j, k, b)),
                        0.0,
                        static_cast<double>(F.y()));
                double z = std::clamp(
                        (static_cast<double>(k)-dt*grid.getW(i, j, k, b)),
                        0.0,
                        static_cast<double>(F.z()));

                std::uint16_t i0 =
                    static_cast<std::uint16_t>(x);
                std::uint16_t i1 =
                    std::clamp(i0 + 1, 1, static_cast<int>(F.x()-1));
                std::uint16_t j0 =
                    static_cast<std::uint16_t>(y);
                std::uint16_t j1 =
                    std::clamp(j0 + 1, 1, static_cast<int>(F.y()-1));
                std::uint16_t k0 =
                    static_cast<std::uint16_t>(z);
                std::uint16_t k1 =
                    std::clamp(k0 + 1, 1, static_cast<int>(F.z()-1));

                const double top = 
                    std::max({F(i0, j0, k0), F(i0, j0, k1), F(i0, j1, k0),
                            F(i0, j1, k1), F(i1, j0, k0), F(i1, j0, k1),
                            F(i1, j1, k0), F(i1, j1, k1)});
                const double bot =
                    std::min({F(i0, j0, k0), F(i0, j0, k1), F(i0, j1, k0),
                            F(i0, j1, k1), F(i1, j0, k0), F(i1, j0, k1),
                            F(i1, j1, k0), F(i1, j1, k1)});

                // Forward step after backward to get error
                x = std::clamp(
                        (static_cast<double>(i)+dt*grid.getU(i, j, k, b)),
                        0.0,
                        static_cast<double>(F.x()));
                y = std::clamp(
                        (static_cast<double>(j)+dt*grid.getV(i, j, k, b)),
                        0.0,
                        static_cast<double>(F.y()));
                z = std::clamp(
                        (static_cast<double>(k)+dt*grid.getW(i, j, k, b)),
                        0.0,
                        static_cast<double>(F.z()));

                const double back = interp(F, x, y, z);
                F(i, j, k) =
                    std::clamp(
                            F(i, j, k) + 0.5 * (Fprev(i, j, k) - back),
                            bot,
                            top
                    );
            }
        }
    }
}

// 3D interpolation in the field F
inline double Advect3D::interp(
        const Field<double, std::uint16_t>& F,
        const double x,
        const double y,
        const double z
    ) const
{
    const std::uint16_t i0 = static_cast<std::uint16_t>(x);
    const std::uint16_t i1 = std::clamp(i0 + 1, 1, static_cast<int>(F.x()-1));
    const std::uint16_t j0 = static_cast<std::uint16_t>(y);
    const std::uint16_t j1 = std::clamp(j0 + 1, 1, static_cast<int>(F.y()-1));
    const std::uint16_t k0 = static_cast<std::uint16_t>(z);
    const std::uint16_t k1 = std::clamp(k0 + 1, 1, static_cast<int>(F.z()-1));

    const double s1 = x - i0;
    const double s0 = 1.0 - s1;
    const double t1 = y - j0;
    const double t0 = 1.0 - t1;
    const double u1 = z - k0;
    const double u0 = 1.0 - u1;

    return s0 * (     t0 * ( u0 * F(i0, j0, k0) + u1 * F(i0, j0, k1) )
                    + t1 * ( u0 * F(i0, j1, k0) + u1 * F(i0, j1, k1) )
                )
         + s1 * (     t0 * ( u0 * F(i1, j0, k0) + u1 * F(i1, j0, k1) )
                    + t1 * ( u0 * F(i1, j1, k0) + u1 * F(i1, j1, k1) )
                );
}

// 2D semi-lagrangian advection, going backward in time to get new values
void Advect2D::advect(
        const StaggeredGrid<double, std::uint16_t>& grid,
        Field<double, std::uint16_t>& F,
        Field<double, std::uint16_t>& Fprev,
        const std::uint8_t b
    )
{
    Fprev = F;
    const double dt = Config::dt * Config::N;
    #pragma omp parallel for
    for (std::uint64_t n = 0; n < grid._surface.maxIt(); ++n)
    {
        const std::uint64_t xy =
            static_cast<std::uint64_t>(grid._surface.x()*grid._surface.y());
        const std::uint64_t m = n % xy;
        const std::uint16_t i = m % grid._surface.x();
        const std::uint16_t j = m / grid._surface.x();
        if (!(F.label(i, j, 0) & SOLID))
        {
            const double x =
                std::clamp((static_cast<double>(i)-dt*grid.getU(i, j, 0, b)),
                            0.0,
                            static_cast<double>(F.x()));
            const double y =
                std::clamp((static_cast<double>(j)-dt*grid.getV(i, j, 0, b)),
                            0.0,
                            static_cast<double>(F.y()));
            F(i, j, 0) = interp(Fprev, x, y);
        }
    }
    if (Config::advection == MACCORMACK)
    {
        // Reverse advection to calculate errors made,
        // than correct the first advection to reduce the errors
        #pragma omp parallel for
        for (std::uint64_t n = 0; n < grid._surface.maxIt(); ++n)
        {
            const std::uint64_t xy =
                static_cast<std::uint64_t>(grid._surface.x()*grid._surface.y());
            const std::uint64_t m = n % xy;
            const std::uint16_t i = m % grid._surface.x();
            const std::uint16_t j = m / grid._surface.x();
            if (!(F.label(i, j, 0) & SOLID))
            {
                double x = std::clamp(
                        (static_cast<double>(i)-dt*grid.getU(i, j, 0, b)),
                        0.0,
                        static_cast<double>(F.x()));
                double y = std::clamp(
                        (static_cast<double>(j)-dt*grid.getV(i, j, 0, b)),
                        0.0,
                        static_cast<double>(F.y()));

                std::uint16_t i0 =
                    static_cast<std::uint16_t>(x);
                std::uint16_t i1 =
                    std::clamp(i0 + 1, 1, static_cast<int>(F.x()-1));
                std::uint16_t j0 =
                    static_cast<std::uint16_t>(y);
                std::uint16_t j1 =
                    std::clamp(j0 + 1, 1, static_cast<int>(F.y()-1));

                const double top = 
                    std::max({F(i0, j0, 0), F(i0, j0, 0), F(i0, j1, 0),
                            F(i0, j1, 0), F(i1, j0, 0), F(i1, j0, 0),
                            F(i1, j1, 0), F(i1, j1, 0)});
                const double bot =
                    std::min({F(i0, j0, 0), F(i0, j0, 0), F(i0, j1, 0),
                            F(i0, j1, 0), F(i1, j0, 0), F(i1, j0, 0),
                            F(i1, j1, 0), F(i1, j1, 0)});

                // Forward step after backward to get error
                x = std::clamp(
                        (static_cast<double>(i)+dt*grid.getU(i, j, 0, b)),
                        0.0,
                        static_cast<double>(F.x()));
                y = std::clamp(
                        (static_cast<double>(j)+dt*grid.getV(i, j, 0, b)),
                        0.0,
                        static_cast<double>(F.y()));

                const double back = interp(F, x, y);
                F(i, j, 0) =
                    std::clamp(
                            F(i, j, 0) + 0.5 * (Fprev(i, j, 0) - back),
                            bot,
                            top
                    );
            }
        }
    }
}

// 2D interpolation in the field F
inline double Advect2D::interp(
        const Field<double, std::uint16_t>& F,
        const double x,
        const double y
    ) const
{
    const std::uint16_t i0 = static_cast<std::uint16_t>(x);
    const std::uint16_t i1 = i0 + 1;
    const std::uint16_t j0 = static_cast<std::uint16_t>(y);
    const std::uint16_t j1 = j0 + 1;

    const double s1 = x - i0;
    const double s0 = 1.0 - s1;
    const double t1 = y - j0;
    const double t0 = 1.0 - t1;

    return s0 * (     t0 * F(i0, j0, 0)
                    + t1 * F(i0, j1, 0)
                )
         + s1 * (     t0 * F(i1, j0, 0)
                    + t1 * F(i1, j1, 0)
                );
}
