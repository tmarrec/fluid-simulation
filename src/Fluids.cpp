#include "Fluids.h"

// Initialise the simulation tools
Fluids::Fluids()
{
    switch (Config::dim)
    {
        case 2:
            _texture = std::vector<std::uint8_t>(
                    _grid._surface.x()*_grid._surface.y()*3
                );
            _advection = std::make_unique<Advect2D>();
            _projection = std::make_unique<Project2D>(_grid);
            break;
        case 3:
            _texture = std::vector<std::uint8_t>(
                    _grid._surface.x()*_grid._surface.y()*_grid._surface.z()
                );
            _advection = std::make_unique<Advect3D>();
            _projection = std::make_unique<Project3D>(_grid);
            break;
    }
}

// Update the simulation by one step
void Fluids::update(const std::uint64_t iteration)
{
    _iteration = iteration;
    
    step();
    if (Config::renderFrames)
    {
        switch (Config::dim)
        {
            case 2:
                updateTexture2D();
                break;
            case 3:
                updateTexture3D();
                break;
        }
    }
}

// Simulation step aka Navier-Stokes solving
void Fluids::step()
{
    for (std::uint16_t k = 0; k < _grid._surface.z(); ++k)
    {
        for (std::uint16_t j = 0; j < _grid._surface.y(); ++j)
        {
            for (std::uint16_t i = 0; i < _grid._surface.x(); ++i)
            {
                double dist = std::sqrt(
                        std::pow(i-_grid._surface.x()/2, 2)+
                        std::pow(j-_grid._surface.y()/3, 2)
                    );
                if (dist < 6 && (k == 1 || k == 2))
                {
                    _grid._surface(i, j, k) = -10;
                    _grid._W(i, j, k) = 5000;
                    _grid._W(i, j, k+1) = 5000;
                    _grid._V(i, j, k) = 7000;
                    _grid._V(i, j, k+1) = 7000;
                }
                else if (dist < 6 &&
                        (k == _grid._surface.z()-2 || k == _grid._surface.z()-3)
                    )
                {
                    _grid._surface(i, j, k) = -10;
                    _grid._W(i, j, k+1) = -5100;
                    _grid._W(i, j, k) = -5100;
                    _grid._V(i, j, k) = 7000;
                    _grid._V(i, j, k-1) = 7000;
                }
                else
                {
                    if (_iteration == 0)
                    {
                        _grid._surface(i, j, k) = 10;
                    }
                }
            }
        }
    }

    // Set labels to fields (inside/outside/..)
    _grid._surface.setLabels(_grid._U, _grid._V, _grid._W);

    // Extrapolate the velocity field
    extrapolate(_grid._U, _grid._UPrev);
    extrapolate(_grid._V, _grid._VPrev);
    extrapolate(_grid._W, _grid._WPrev);

    // Advect level-set everywhere using the fully extrapolated velocity
    _advection->advect(_grid, _grid._surface, _grid._surfacePrev, 0);
    redistancing(8, _grid._surface, _grid._surfacePrev);

    // Advect velocity everywhere using the fully extrapolated velocity
    _advection->advect(_grid, _grid._U, _grid._UPrev, 1);
    _advection->advect(_grid, _grid._V, _grid._VPrev, 2);
    _advection->advect(_grid, _grid._W, _grid._WPrev, 3);

    // Set labels to fields (inside/outside/..)
    _grid._surface.setLabels(_grid._U, _grid._V, _grid._W);

    // Add external forces
    addForces();

    // Tag the cells that are inside the liquids and assign integer labels
    _grid.tagActiveCells();

    // Ensure incompressibility
    _projection->project();
}

// Add forces to the velocity field (gravity for example)
void Fluids::addForces()
{
    const double G = 150.0;
    for (std::uint16_t k = 0; k < _grid._V.z(); ++k)
    {
        for (std::uint16_t j = 0; j < _grid._V.y(); ++j)
        {
            for (std::uint16_t i = 0; i < _grid._V.x(); ++i)
            {
                if (!(_grid._V.label(i, j, k) & SOLID))
                {
                    _grid._V(i, j, k) -= G;
                }
            }
        }
    }
}

// Extrapolate the field F, nbIte times
void Fluids::extrapolate(
        Field<double, std::uint16_t>& F,
        Field<double, std::uint16_t>& Ftemp,
        std::uint16_t nbIte
    ) const
{
    // Use Ftemp to store each extrapolation step
    // For each cell, if its neighbors are either LIQUID, SOLID or EXTRAPOLATED
    // then set it to EXTRAPOLATED with the average value of its valid neighbors
    std::uint16_t it = 0;
    std::uint64_t nbNeg = 0;
    do
    {
        nbNeg = 0;
        #pragma omp parallel for
        for (std::uint64_t n = 0; n < F.maxIt(); ++n)
        {
            const std::uint64_t xy = static_cast<std::uint64_t>(F.x()*F.y());
            const std::uint64_t m = n % xy;
            const std::uint16_t i = m % F.x();
            const std::uint16_t j = m / F.x();
            const std::uint16_t k = n / xy;

            if (F.label(i, j, k) & LIQUID)
            {
                Ftemp(i, j, k) = F(i, j, k);
                Ftemp.label(i, j, k) = LIQUID;
            }
            else if (F.label(i, j, k) & SOLID)
            {
                Ftemp(i, j, k) = F(i, j, k);
                Ftemp.label(i, j, k) = SOLID;
            }
            else if (F.label(i, j, k) & EXTRAPOLATED)
            {
                Ftemp(i, j, k) = F(i, j, k);
                Ftemp.label(i, j, k) = Ftemp.label(i, j, k) | EXTRAPOLATED;
            }
            else
            {
                std::uint8_t nbNeighbors = 0;
                double value = 0.0;
                if (i < _grid._surface.x()-1 && F.checked(i+1, j, k))
                {
                    nbNeighbors++;
                    value += F(i+1, j, k);
                }
                if (i > 0 && F.checked(i-1, j, k))
                {
                    nbNeighbors++;
                    value += F(i-1, j, k);
                }
                if (j < _grid._surface.y()-1 && F.checked(i, j+1, k))
                {
                    nbNeighbors++;
                    value += F(i, j+1, k);
                }
                if (j > 0 && F.checked(i, j-1, k))
                {
                    nbNeighbors++;
                    value += F(i, j-1, k);
                }
                if (k < _grid._surface.z()-1 && F.checked(i, j, k+1))
                {
                    nbNeighbors++;
                    value += F(i, j, k+1);
                }
                if (k > 0 && F.checked(i, j, k-1))
                {
                    nbNeighbors++;
                    value += F(i, j, k-1);
                }
                if (nbNeighbors > 0)
                {
                    nbNeg++;
                    Ftemp(i, j, k) = value/nbNeighbors;
                    Ftemp.label(i, j, k) = Ftemp.label(i, j, k) | EXTRAPOLATED;
                }
            }
        }
        F = Ftemp;
        if (nbIte > 0 && ++it == nbIte)
        {
            return;
        }
    } while (nbNeg > 0);
}

// Try to force the gradient norm of the level-set to be equal to 1
void Fluids::redistancing(
        const std::uint64_t nbIte,
        Field<double, std::uint16_t>& field,
        Field<double, std::uint16_t>& fieldTemp
    ) const
{
    const double dx = 1.0/Config::N;
    auto F = field;
    auto QNew = field;
    for (std::int16_t k = 0; k < field.z(); ++k)
    {
        for (std::int16_t j = 0; j < field.y(); ++j)
        {
            for (std::int16_t i = 0; i < field.x(); ++i)
            {
                if ((i+1 < field.x() &&
                        !((field(i, j, k) >= 0) ^ (field(i+1, j, k) < 0))) ||
                    (i-1 >= 0 &&
                        !((field(i, j, k) >= 0) ^ (field(i-1, j, k) < 0))) ||
                    (j+1 < field.y() &&
                        !((field(i, j, k) >= 0) ^ (field(i, j+1, k) < 0))) ||
                    (j-1 >= 0 &&
                        !((field(i, j, k) >= 0) ^ (field(i, j-1, k) < 0))) ||
                    (k+1 < field.z() &&
                        !((field(i, j, k) >= 0) ^ (field(i, j, k+1) < 0))) ||
                    (k-1 >= 0 &&
                        !((field(i, j, k) >= 0) ^ (field(i, j, k-1) < 0))))
                {
                    F(i, j, k) = 1;
                    F.label(i, j, k) = LIQUID;
                }
                else
                {
                    F(i, j, k) = 0;
                    F.label(i, j, k) = EMPTY;
                }
            }
        }
    }

    extrapolate(F, fieldTemp, 2);

    std::uint8_t dist = 2;
    for (std::uint16_t k = 0; k < field.z(); ++k)
    {
        for (std::uint16_t j = 0; j < field.y(); ++j)
        {
            for (std::uint16_t i = 0; i < field.x(); ++i)
            {
                if (F(i, j, k) == 0)
                {
                    QNew(i, j, k) = dist * (field(i, j, k) <= 0 ? -1 : 1);
                }
            }
        }
    }
    for (std::uint16_t k = 0; k < field.z(); ++k)
    {
        for (std::uint16_t j = 0; j < field.y(); ++j)
        {
            for (std::uint16_t i = 0; i < field.x(); ++i)
            {
                if (std::abs(QNew(i, j, k)) > dist)
                {
                    QNew(i, j, k) = dist * (field(i, j, k) <= 0 ? -1 : 1);
                }
            }
        }
    }
    field = QNew;
    Field<double, std::uint16_t> n = field;
    Field<double, std::uint16_t> Ssf = field;
    // Init smoothing function
    for (std::uint16_t k = 0; k < field.z()-0; ++k)
    {
        for (std::uint16_t j = 0; j < field.y()-0; ++j)
        {
            for (std::uint16_t i = 0; i < field.x()-0; ++i)
            {
                const double O0 = field(i, j, k);
                Ssf(i, j, k) =
                    O0 / (std::sqrt(std::pow(O0, 2) + std::pow(0.5, 2)));
            }
        }
    }

    // Step forward in fictious time
    for (std::uint64_t relaxit = 0; relaxit < nbIte; ++relaxit)
    {
        for (std::uint16_t k = 0; k < field.z(); ++k)
        {
            for (std::uint16_t j = 0; j < field.y(); ++j)
            {
                for (std::uint16_t i = 0; i < field.x(); ++i)
                {
                    if (F(i, j, k) == 1)
                    {
                        const double gO = field.gradLength(i, j, k);
                        n(i, j, k) =
                            (0.5 * dx * (- Ssf(i, j, k) * (gO - 1.0))) +
                            field(i, j, k);
                    }
                }
            }
        }
        field = n;
    }
}

// Export the level-set to a 3D texture
void Fluids::updateTexture3D()
{
    std::uint64_t it = 0;
    for (std::uint16_t k = 0; k < _grid._surface.z(); ++k)
    {
        for (std::uint16_t j = 0; j < _grid._surface.y(); ++j)
        {
            for (std::uint16_t i = 0; i < _grid._surface.x(); ++i)
            {
                _texture[it] = _grid._surface(i, j, k) < 0 ? 32 : 0;
                it++;
            }
        }
    }
}

// Export the level-set to a 2D texture
void Fluids::updateTexture2D()
{
    std::uint64_t it = 0;
    for (std::uint16_t i = 0; i < _grid._surface.x(); ++i)
    {
        for (std::uint16_t j = 0; j < _grid._surface.y(); ++j)
        {
            const double implicit = _grid._surface(i, j, 0) <= 0 ? -255 : 255;
            double p = 0;
            double pneg = 0;
            if (_grid._pressure(i, j, 0) < 0)
            {
                pneg = _grid._pressure(i, j, 0)*(-1);
            }
            else
            {
                p = _grid._pressure(i, j, 0);
            }
            if (implicit <= 0)
            {
                _texture[it*3+0] = std::clamp(p*150, 0.0, 255.0);
                _texture[it*3+1] = std::clamp(pneg*150, 0.0, 255.0);
                _texture[it*3+2] = std::clamp(-implicit*50, 0.0, 255.0);
            }
            else
            {
                _texture[it*3+0] = 0;
                _texture[it*3+1] = 0;
                _texture[it*3+2] = 0;
            }
            it++;
        }
    }
}

// Used to render velocity field in 2D
const std::vector<double>& Fluids::X() const
{
    return _grid._U.data();
}

// Used to render velocity field in 2D
const std::vector<double>& Fluids::Y() const
{
    return _grid._V.data();
}

// Used to render velocity field in 2D
// Return if the cell (i,j,k) is active or not
bool Fluids::isCellActive(
        const std::uint16_t i,
        const std::uint16_t j,
        const std::uint16_t k
    ) const
{
    return _grid._pressureID(i, j, k) > 0;
}

// Used to render velocity field in 2D
const Field<double, std::uint16_t>& Fluids::surface() const
{
    return _grid._surface;
}

// Used to render velocity field in 2D
const std::vector<std::uint8_t>& Fluids::texture() const
{
    return _texture;
}


