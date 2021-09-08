#include "Fluids.h"
#include "ConjugateGradient.h"
#include <iomanip>

Fluids::Fluids()
{
    if (Config::dim == 2)
    {
        _texture = std::vector<std::uint8_t>(_grid._surface.x()*_grid._surface.y()*3);
    }
    else if (Config::dim == 3)
    {
        _texture = std::vector<std::uint8_t>(_grid._surface.x()*_grid._surface.y()*_grid._surface.z());
    }
}

void Fluids::update([[maybe_unused]] std::uint64_t iteration)
{
    _iteration = iteration;
    
    step();

    updateTexture();

    //std::cin.get();
}

void Fluids::step()
{
    _grid._U.resetPos();
    _grid._V.resetPos();
    _grid._W.resetPos();

    for (std::uint16_t k = 0; k < _grid._surface.z(); ++k)
    {
        for (std::uint16_t j = 0; j < _grid._surface.y(); ++j)
        {
            for (std::uint16_t i = 0; i < _grid._surface.x(); ++i)
            {
                if (_grid._surface(i,j,k) < 0)
                {
                    _grid._U.label(i,j,k) = LIQUID;
                    _grid._U.label(i+1,j,k) = LIQUID;
                    _grid._V.label(i,j,k) = LIQUID;
                    _grid._V.label(i,j+1,k) = LIQUID;
                    _grid._W.label(i,j,k) = LIQUID;
                    _grid._W.label(i,j,k+1) = LIQUID;
                }
            }
        }
    }

    for (std::uint16_t k = 0; k < std::max({_grid._U.z(), _grid._V.z(), _grid._W.z()}); ++k)
    {
        for (std::uint16_t j = 0; j < std::max({_grid._U.y(), _grid._V.y(), _grid._W.y()}); ++j)
        {
            for (std::uint16_t i = 0; i < std::max({_grid._U.x(), _grid._V.x(), _grid._W.x()}); ++i)
            {
                if (k < _grid._U.z() && j < _grid._U.y() && i < _grid._U.x())
                {
                    if (i == 0 || i == _grid._U.x()-1 || j == 0 || j == _grid._U.y()-1 || k == 0 || k == _grid._U.z()-1)
                    {
                        _grid._U.label(i,j,k) = SOLID;
                    }
                }
                if (k < _grid._V.z() && j < _grid._V.y() && i < _grid._V.x())
                {
                    if (i == 0 || i == _grid._V.x()-1 || j == 0 || j == _grid._V.y()-1 || k == 0 || k == _grid._V.z()-1)
                    {
                        _grid._V.label(i,j,k) = SOLID;
                    }
                }
                if (k < _grid._W.z() && j < _grid._W.y() && i < _grid._W.x())
                {
                    if (i == 0 || i == _grid._W.x()-1 || j == 0 || j == _grid._W.y()-1 || k == 0 || k == _grid._W.z()-1)
                    {
                        _grid._W.label(i,j,k) = SOLID;
                    }
                }
            }
        }
    }

    // Extrapolate the velocity field
    extrapolate(_grid._U, _grid._UPrev);
    extrapolate(_grid._V, _grid._VPrev);
    extrapolate(_grid._W, _grid._WPrev);

    // Advect level-set everywhere using the fully extrapolated velocity
    advect(_grid._surface, _grid._surfacePrev, 0);
    //redistancing(4, _grid._surface, _grid._surfacePrev);

    // Advect velocity everywhere using the fully extrapolated velocity
    advect(_grid._U, _grid._UPrev, 1);
    advect(_grid._V, _grid._VPrev, 2);
    advect(_grid._W, _grid._WPrev, 3);

    // Add external forces
    addForces();

    // Tag the cells that are inside the liquids and assign integer labels
    _grid.tagActiveCells();

    // Ensure incompressibility
    project();
}

void Fluids::addForces()
{
    {
        for (std::uint16_t k = 0; k < _grid._surface.z(); ++k)
        {
            for (std::uint16_t j = 0; j < _grid._surface.y(); ++j)
            {
                for (std::uint16_t i = 0; i < _grid._surface.x(); ++i)
                {
                    double dist = std::sqrt(std::pow(i-_grid._surface.x()/2,2)+std::pow(j-_grid._surface.y()/2,2)+std::pow(k-_grid._surface.z()/2,2));
                    if (dist < 8)
                    {
                        if (_iteration < 256)
                        _grid._surface(i,j,k) = -dist;
                    }
                    else
                    {
                        if (_iteration == 0)
                        _grid._surface(i,j,k) = 10;
                    }
                    //_grid._surface(i,j,k) = 10;
                }
            }
        }
    }
    const double G = -15.0;
    std::uint64_t s = 8;
    for (std::uint16_t k = 0; k < _grid._V.z(); ++k)
    {
        for (std::uint16_t j = 0; j < _grid._V.y(); ++j)
        {
            for (std::uint16_t i = 0; i < _grid._V.x(); ++i)
            {
                if (!(_grid._V.label(i,j,k) & SOLID))
                {
                    _grid._V(i,j,k) += G;
                    
                    /*
                    std::uint16_t cubeSize = 4;
                    if (
                        ((i >= 5 && i < 5+cubeSize) && (j >= 24 && j < 24+cubeSize) && (k >= 5 && k < 5+cubeSize))
                    )
                    {
                        _grid._surface(i,j,k) = -10;
                        _grid._V(i,j,k) = -100;
                        _grid._V(i,j+1,k) = -100;
                    }
                    */
                }
            }
        }
    }
}

inline double Fluids::interp(const Field<double,std::uint16_t>& F, double x, double y, double z) const
{
    const std::uint16_t i0 = static_cast<std::uint16_t>(x);
    const std::uint16_t i1 = i0 + 1;
    const std::uint16_t j0 = static_cast<std::uint16_t>(y);
    const std::uint16_t j1 = j0 + 1;
    const std::uint16_t k0 = static_cast<std::uint16_t>(z);
    const std::uint16_t k1 = k0 + 1;

    const double s1 = x - i0;
    const double s0 = 1.0 - s1;
    const double t1 = y - j0;
    const double t0 = 1.0 - t1;
    const double u1 = z - k0;
    const double u0 = 1.0 - u1;

    return s0 * (     t0 * ( u0 * F(i0,j0,k0) + u1 * F(i0,j0,k1) )
                    + t1 * ( u0 * F(i0,j1,k0) + u1 * F(i0,j1,k1) )
                )
         + s1 * (     t0 * ( u0 * F(i1,j0,k0) + u1 * F(i1,j0,k1) )
                    + t1 * ( u0 * F(i1,j1,k0) + u1 * F(i1,j1,k1) )
                );
}

void Fluids::extrapolate(Field<double,std::uint16_t>& F, Field<double,std::uint16_t>& Ftemp, std::uint16_t nbIte)
{
    std::uint16_t it = 0;
    std::uint64_t nbNeg = 0;
    do
    {
        nbNeg = 0;
        #pragma omp parallel for
        for (std::uint64_t n = 0; n < F.maxIt(); ++n)
        {
            const std::uint64_t xy = F.x()*F.y();
            const std::uint64_t m = n % xy;
            const std::uint16_t i = m % F.x();
            const std::uint16_t j = m / F.x();
            const std::uint16_t k = n / xy;

            if (F.label(i,j,k) & LIQUID)
            {
                Ftemp(i,j,k) = F(i,j,k);
                Ftemp.label(i,j,k) = LIQUID;
            }
            else if (F.label(i,j,k) & SOLID)
            {
                Ftemp(i,j,k) = F(i,j,k);
                Ftemp.label(i,j,k) = SOLID;
            }
            else if (F.label(i,j,k) & EXTRAPOLATED)
            {
                Ftemp(i,j,k) = F(i,j,k);
                Ftemp.label(i,j,k) = Ftemp.label(i,j,k) | EXTRAPOLATED;
            }
            else
            {
                std::uint8_t nbNeighbors = 0;
                double value = 0.0;
                if (i < Config::N-1 && F.checked(i+1,j,k))
                {
                    nbNeighbors++;
                    value += F(i+1,j,k);
                }
                if (i > 0 && F.checked(i-1,j,k))
                {
                    nbNeighbors++;
                    value += F(i-1,j,k);
                }
                if (j < Config::N-1 && F.checked(i,j+1,k))
                {
                    nbNeighbors++;
                    value += F(i,j+1,k);
                }
                if (j > 0 && F.checked(i,j-1,k))
                {
                    nbNeighbors++;
                    value += F(i,j-1,k);
                }
                if (k < Config::N-1 && F.checked(i,j,k+1))
                {
                    nbNeighbors++;
                    value += F(i,j,k+1);
                }
                if (k > 0 && F.checked(i,j,k-1))
                {
                    nbNeighbors++;
                    value += F(i,j,k-1);
                }
                if (nbNeighbors > 0)
                {
                    nbNeg++;
                    Ftemp(i,j,k) = value/nbNeighbors;
                    Ftemp.label(i,j,k) = Ftemp.label(i,j,k) | EXTRAPOLATED;
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

void Fluids::redistancing(const std::uint64_t nbIte, Field<double, std::uint16_t>& field, Field<double, std::uint16_t>& fieldTemp)
{
    WARNING("Redistancing probably wrong (0.5 or 0.333???)");
    const double dx = 1.0/Config::N;
    auto F = field;
    auto QNew = field;
    for (std::uint16_t k = 0; k < field.z(); ++k)
    {
        for (std::uint16_t j = 0; j < field.y(); ++j)
        {
            for (std::uint16_t i = 0; i < field.x(); ++i)
            {
                if (    (i+1 < field.x() && !((field(i,j,k) >= 0) ^ (field(i+1,j,k) < 0))) ||
                        (i-1 >= 0 && !((field(i,j,k) >= 0) ^ (field(i-1,j,k) < 0))) ||
                        (j+1 < field.y() && !((field(i,j,k) >= 0) ^ (field(i,j+1,k) < 0))) ||
                        (j-1 >= 0 && !((field(i,j,k) >= 0) ^ (field(i,j-1,k) < 0))) ||
                        (k+1 < field.z() && !((field(i,j,k) >= 0) ^ (field(i,j,k+1) < 0))) ||
                        (k-1 >= 0 && !((field(i,j,k) >= 0) ^ (field(i,j,k-1) < 0)))   )
                {
                    F(i,j,k) = 1;
                    F.label(i,j,k) = LIQUID;
                }
                else
                {
                    F(i,j,k) = 0;
                    F.label(i,j,k) = EMPTY;
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
                if (F(i,j,k) == 0)
                {
                    QNew(i,j,k) = dist * (field(i,j,k) <= 0 ? -1 : 1);
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
                if (std::abs(QNew(i,j,k)) > dist)
                {
                    QNew(i,j,k) = dist * (field(i,j,k) <= 0 ? -1 : 1);
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
                const double O0 = field(i,j,k);
                Ssf(i,j,k) = O0 / (std::sqrt(std::pow(O0, 2) + std::pow(0.5, 2)));
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
                    if (F(i,j,k) == 1)
                    {
                        const double gO = field.gradLength(i,j,k);
                        n(i,j,k) = (0.5 * dx * (- Ssf(i,j,k) * (gO - 1.0))) + field(i,j,k);
                    }
                }
            }
        }
        field = n;
    }
}

void Fluids::diffuse(Field<double,std::uint16_t>& F, const Field<double,std::uint16_t>& Fprev, const std::uint8_t b, const Laplacian& A)
{
    /*
    Eigen::VectorXd x(A.diag.size());
    ConjugateGradient(A, x, Fprev.vec(), _solverType);
    F.setFromVec(x);
    */
}


void Fluids::advect(Field<double,std::uint16_t>& F, Field<double,std::uint16_t>& Fprev, const std::uint8_t b) const
{
    Fprev = F;
	const double dt = Config::dt * Config::N;
    for (std::uint16_t k = 0; k < F.z(); ++k)
    {
        for (std::uint16_t j = 0; j < F.y(); ++j)
        {
            for (std::uint16_t i = 0; i < F.x(); ++i)
            {
                if (!(F.label(i,j,k) & SOLID))
                {
                    const double posx = static_cast<double>(i);
                    const double posy = static_cast<double>(j);
                    const double posz = static_cast<double>(k);

                    const double x = std::clamp((posx-dt*_grid.getU(i,j,k,b)), 0.0, static_cast<double>(F.x()+0.0));
                    const double y = std::clamp((posy-dt*_grid.getV(i,j,k,b)), 0.0, static_cast<double>(F.y()+0.0));
                    const double z = std::clamp((posz-dt*_grid.getW(i,j,k,b)), 0.0, static_cast<double>(F.z()+0.0));

                    F(i,j,k) = interp(Fprev, x, y, z);
                }
            }
        }
    }
}

void Fluids::preparePressureSolving(Laplacian& laplacian, Eigen::VectorXd& b)
{
    laplacian.A = Eigen::SparseMatrix<double>(_grid.activeCellsNb(), _grid.activeCellsNb());
    laplacian.A.reserve(Eigen::VectorXi::Constant(_grid.activeCellsNb(), 7));
    _grid._Adiag.reset();
    _grid._Ax.reset();
    _grid._Ay.reset();
    _grid._Az.reset();
    _grid._precon.reset();

    const double scale = 1;
    for (std::uint64_t n = 0; n < _grid._surface.maxIt(); ++n)
    {
        const std::uint64_t xy = _grid._surface.x()*_grid._surface.y();
        const std::uint64_t m = n % xy;
        const std::uint16_t i = m % _grid._surface.x();
        const std::uint16_t j = m / _grid._surface.x();
        const std::uint16_t k = n / xy;

        std::uint64_t id = _grid._pressureID(i,j,k);
        if (id > 0)
        {
            const std::uint64_t realID = id-1;
            std::uint64_t neibID = 0;
            if (i > 0 && (neibID = _grid._pressureID(i-1,j,k)) > 0)
            {
                laplacian.A.coeffRef(realID, neibID-1) = -scale;
                _grid._Adiag(i,j,k) += scale;
            }
            if (i+1 < Config::N && (neibID = _grid._pressureID(i+1,j,k)) > 0)
            {
                laplacian.A.coeffRef(realID, neibID-1) = -scale;
                _grid._Adiag(i,j,k) += scale;
                _grid._Ax(i,j,k) = -scale;
            }
            else if (i+1 < Config::N && _grid._surface.label(i+1,j,k) & EMPTY)
            {
                _grid._Adiag(i,j,k) += scale;
            }
            if (j > 0 && (neibID = _grid._pressureID(i,j-1,k)) > 0)
            {
                laplacian.A.coeffRef(realID, neibID-1) = -scale;
                _grid._Adiag(i,j,k) += scale;
            }
            if (j+1 < Config::N && (neibID = _grid._pressureID(i,j+1,k)) > 0)
            {
                laplacian.A.coeffRef(realID, neibID-1) = -scale;
                _grid._Adiag(i,j,k) += scale;
                _grid._Ay(i,j,k) = -scale;
            }
            else if (j+1 < Config::N && _grid._surface.label(i,j+1,k) & EMPTY)
            {
                _grid._Adiag(i,j,k) += scale;
            }
            if (k > 0 && (neibID = _grid._pressureID(i,j,k-1)) > 0)
            {
                laplacian.A.coeffRef(realID, neibID-1) = -scale;
                _grid._Adiag(i,j,k) += scale;
            }
            if (k+1 < Config::N && (neibID = _grid._pressureID(i,j,k+1)) > 0)
            {
                laplacian.A.coeffRef(realID, neibID-1) = -scale;
                _grid._Adiag(i,j,k) += scale;
                _grid._Az(i,j,k) = -scale;
            }
            else if (k+1 < Config::N && _grid._surface.label(i,j,k+1) & EMPTY)
            {
                _grid._Adiag(i,j,k) += scale;
            }

            std::uint8_t nonSolidNeib = 6;
            if (i == 0)
                nonSolidNeib--;
            if (j == 0)
                nonSolidNeib--;
            if (k == 0)
                nonSolidNeib--;
            if (i == _grid._surface.x()-1)
                nonSolidNeib--;
            if (j == _grid._surface.y()-1)
                nonSolidNeib--;
            if (k == _grid._surface.z()-1)
                nonSolidNeib--;
            laplacian.A.coeffRef(realID, realID) = nonSolidNeib;
            b.coeffRef(realID) = div(i,j,k);
        }
    }
}

inline double Fluids::div(const std::uint16_t i, const std::uint16_t j, const std::uint16_t k) const
{
    const double h = 1.0/Config::N;
    double Udiv = _grid._U(i+1,j,k) - _grid._U(i,j,k);
    double Vdiv = _grid._V(i,j+1,k) - _grid._V(i,j,k);
    double Zdiv = _grid._V(i,j,k+1) - _grid._V(i,j,k);
    double rhs = - h * (Udiv + Vdiv + Zdiv);
    return rhs;
}

inline double Fluids::pressureAt(const std::uint16_t i, const std::uint16_t j, const std::uint16_t k, const Eigen::VectorXd x, const std::uint64_t l) const
{
    return _grid._pressure(i,j,k);
}

void Fluids::project()
{
    if (_grid.activeCellsNb() > 0)
    {
        Eigen::VectorXd x (_grid.activeCellsNb());
        Eigen::VectorXd b (_grid.activeCellsNb());

        Laplacian A {};
        preparePressureSolving(A, b);
        ConjugateGradient(A, x, b, _grid);

        _grid._pressure.reset();

        std::uint64_t id;
        for (std::uint16_t k = 0; k < std::max({_grid._U.z(), _grid._V.z(), _grid._W.z()}); ++k)
        {
            for (std::uint16_t j = 0; j < std::max({_grid._U.y(), _grid._V.y(), _grid._W.y()}); ++j)
            {
                for (std::uint16_t i = 0; i < std::max({_grid._U.x(), _grid._V.x(), _grid._W.x()}); ++i)
                {
                    if (k < _grid._surface.z() && j < _grid._surface.y() && i < _grid._surface.x() && (id = _grid._pressureID(i,j,k)) > 0)
                    {
                        _grid._pressure(i,j,k) = x(id-1);
                    }
                    if (k < _grid._U.z() && j < _grid._U.y() && i < _grid._U.x() && _grid._U.label(i,j,k) == LIQUID)
                    {
                        _grid._U(i,j,k) -= Config::N*(pressureAt(i,j,k,x,0) - pressureAt(i-1,j,k,x,0));
                    }
                    if (k < _grid._V.z() && j < _grid._V.y() && i < _grid._V.x() && _grid._V.label(i,j,k) == LIQUID)
                    {
                        _grid._V(i,j,k) -= Config::N*(pressureAt(i,j,k,x,0) - pressureAt(i,j-1,k,x,0));
                    }
                    if (k < _grid._W.z() && j < _grid._W.y() && i < _grid._W.x() && _grid._W.label(i,j,k) == LIQUID)
                    {
                        _grid._W(i,j,k) -= Config::N*(pressureAt(i,j,k,x,0) - pressureAt(i,j,k-1,x,0));
                    }
                }
            }
        }
    }
}

void Fluids::updateTexture()
{
    std::uint64_t it = 0;
    //double sum = 0;
    for (std::uint16_t k = 0; k < _grid._surface.z(); ++k)
    {
        for (std::uint16_t j = 0; j < _grid._surface.y(); ++j)
        {
            for (std::uint16_t i = 0; i < _grid._surface.x(); ++i)
            {
                _texture[it] = _grid._surface(i,j,k) < 0 ? 32 : 0;
                it++;
            }
        }
	}
}

const std::vector<std::uint8_t>& Fluids::texture() const
{
    return _texture;
}

void Fluids::writeVolumeFile([[maybe_unused]] const std::uint64_t iteration)
{
}

const std::vector<double>& Fluids::X() const
{
    return _grid._U.data();
}

const std::vector<double>& Fluids::Y() const
{
    return _grid._V.data();
}

const std::vector<double>& Fluids::Z() const
{
    return _grid._W.data();
}

bool Fluids::isCellActive(std::uint16_t i, std::uint16_t j, std::uint16_t k) const
{
    return _grid._pressureID(i,j,k) > 0;
}
