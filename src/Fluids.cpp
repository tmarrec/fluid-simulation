#include "Fluids.h"
#include "ConjugateGradient.h"
#include <iomanip>


Fluids::Fluids()
{
    //initLaplacians(Config::N, _dt, _viscosity, _diffusion, _laplacianViscosityX, _laplacianViscosityY, _laplacianProject, _laplacianDiffuse);
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

    for (std::uint16_t j = 0; j < _grid._surface.y()+1; ++j)
    {
        for (std::uint16_t i = 0; i < _grid._surface.x()+1; ++i)
        {
            if (_grid._surface(i,j) < 0)
            {
                _grid._U.label(i,j) = LIQUID;
                _grid._U.label(i+1,j) = LIQUID;
                _grid._V.label(i,j) = LIQUID;
                _grid._V.label(i,j+1) = LIQUID;
            }
        }
    }

    for (std::uint16_t j = 0; j < _grid._surface.y()+1; ++j)
    {
        for (std::uint16_t i = 0; i < _grid._surface.x()+1; ++i)
        {
            if (j < _grid._U.y() && (i == 0 || i == _grid._surface.x()))
            {
                _grid._U.label(i,j) = SOLID;
            }
            if (i < _grid._V.x() && (j == 0 || j == _grid._surface.y()))
            {
                _grid._V.label(i,j) = SOLID;
            }
        }
    }

    // Extrapolate the velocity field
    extrapolate(_grid._U, _grid._UPrev);
    extrapolate(_grid._V, _grid._VPrev);

    // Advect level-set everywhere using the fully extrapolated velocity
    advect(_grid._surface, _grid._surfacePrev, 0);
    redistancing(4, _grid._surface, _grid._surfacePrev);

    // Advect velocity everywhere using the fully extrapolated velocity
    advect(_grid._U, _grid._UPrev, 1);
    advect(_grid._V, _grid._VPrev, 2);

    // Add external forces
    addForces();

    // Tag the cells that are inside the liquids and assign integer labels
    _grid.tagActiveCells();

    // Ensure incompressibility
    project();
}

void Fluids::addForces()
{
    if (_iteration == 0)
    {
        for (std::uint16_t j = 0; j < Config::N; ++j)
        {
            for (std::uint16_t i = 0; i < Config::N; ++i)
            {
                _grid._surface(i,j) = 10;
            }
        }
    }
    const double G = 15.0;
    std::uint64_t s = 64;
    for (std::uint16_t j = 0; j < _grid._V.y(); ++j)
    {
        for (std::uint16_t i = 0; i < _grid._V.x(); ++i)
        {
            if (!(_grid._V.label(i,j) & SOLID))
            {
                _grid._V(i,j) += G;
                std::uint16_t cubeSize = 32;
                if (
                    ((i >= s && i < s+cubeSize) && (j >= 64 && j < 64+cubeSize))
                )
                {
                    _grid._surface(i,j) = -10;
                    _grid._V(i,j) = 500;
                    _grid._V(i,j+1) = 500;
                }
            }
        }
    }
}

inline double Fluids::interp(const Field<double,std::uint16_t>& F, double x, double y) const
{
    const std::uint16_t i0 = static_cast<std::uint16_t>(x);
    const std::uint16_t i1 = i0 + 1;
    const std::uint16_t j0 = static_cast<std::uint16_t>(y);
    const std::uint16_t j1 = j0 + 1;

    const double s1 = x - i0;
    const double s0 = 1.0 - s1;
    const double t1 = y - j0;
    const double t0 = 1.0 - t1;

    const double vA = F(i0,j0);
    const double vB = F(i0,j1);
    const double vC = F(i1,j0);
    const double vD = F(i1,j1);

    return s0*(t0*vA+t1*vB)+s1*(t0*vC+t1*vD);
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
            const std::uint64_t m = n % F.maxIt();
            const std::uint16_t i = m % F.x();
            const std::uint16_t j = m / F.x();

            if (F.label(i,j) & LIQUID)
            {
                Ftemp(i,j) = F(i,j);
                Ftemp.label(i,j) = LIQUID;
            }
            else if (F.label(i,j) & SOLID)
            {
                Ftemp(i,j) = F(i,j);
                Ftemp.label(i,j) = SOLID;
            }
            else if (F.label(i,j) & EXTRAPOLATED)
            {
                Ftemp(i,j) = F(i,j);
                Ftemp.label(i,j) = Ftemp.label(i,j) | EXTRAPOLATED;
            }
            else
            {
                std::uint8_t nbNeighbors = 0;
                double value = 0.0;
                if (i < Config::N-1 && F.checked(i+1,j))
                {
                    nbNeighbors++;
                    value += F(i+1,j);
                }
                if (i > 0 && F.checked(i-1,j))
                {
                    nbNeighbors++;
                    value += F(i-1,j);
                }
                if (j < Config::N-1 && F.checked(i,j+1))
                {
                    nbNeighbors++;
                    value += F(i,j+1);
                }
                if (j > 0 && F.checked(i,j-1))
                {
                    nbNeighbors++;
                    value += F(i,j-1);
                }
                if (nbNeighbors > 0)
                {
                    nbNeg++;
                    Ftemp(i,j) = value/nbNeighbors;
                    Ftemp.label(i,j) = Ftemp.label(i,j) | EXTRAPOLATED;
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
    const double dx = 1.0/Config::N;
    auto F = field;
    auto QNew = field;
    for (std::uint16_t j = 0; j < field.y(); ++j)
    {
        for (std::uint16_t i = 0; i < field.x(); ++i)
        {
            if (    (i+1 < field.x() && !((field(i,j) >= 0) ^ (field(i+1,j) < 0))) ||
                    (i-1 >= 0 && !((field(i,j) >= 0) ^ (field(i-1,j) < 0))) ||
                    (j+1 < field.y() && !((field(i,j) >= 0) ^ (field(i,j+1) < 0))) ||
                    (j-1 >= 0 && !((field(i,j) >= 0) ^ (field(i,j-1) < 0)))   )
            {
                F(i,j) = 1;
                F.label(i,j) = LIQUID;
            }
            else
            {
                F(i,j) = 0;
                F.label(i,j) = EMPTY;
            }
        }
    }
    extrapolate(F, fieldTemp, 2);

    std::uint8_t dist = 2;
    for (std::uint16_t j = 0; j < field.y(); ++j)
    {
        for (std::uint16_t i = 0; i < field.x(); ++i)
        {
            if (F(i,j) == 0)
            {
                QNew(i,j) = dist * (field(i,j) <= 0 ? -1 : 1);
            }
        }
    }
    for (std::uint16_t j = 0; j < field.y(); ++j)
    {
        for (std::uint16_t i = 0; i < field.x(); ++i)
        {
            if (std::abs(QNew(i,j)) > dist)
            {
                QNew(i,j) = dist * (field(i,j) <= 0 ? -1 : 1);
            }
        }
    }
    field = QNew;
    Field<double, std::uint16_t> n = field;
    Field<double, std::uint16_t> Ssf = field;
    // Init smoothing function
    for (std::uint16_t j = 0; j < field.y()-0; ++j)
    {
        for (std::uint16_t i = 0; i < field.x()-0; ++i)
        {
            const double O0 = field(i,j);
            Ssf(i,j) = O0 / (std::sqrt(std::pow(O0, 2) + std::pow(0.5, 2)));
        }
    }

    // Step forward in fictious time
    for (std::uint64_t relaxit = 0; relaxit < nbIte; ++relaxit)
    {
        for (std::uint16_t j = 0; j < field.y(); ++j)
        {
            for (std::uint16_t i = 0; i < field.x(); ++i)
            {
                if (F(i,j) == 1)
                {
                    const double gO = field.gradLength(i, j);
                    n(i,j) = (0.5 * dx * (- Ssf(i,j) * (gO - 1.0))) + field(i,j);
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
    for (std::uint16_t j = 0; j < F.y(); ++j)
    {
        for (std::uint16_t i = 0; i < F.x(); ++i)
        {
            if (!(F.label(i,j) & SOLID))
            {
                const double posx = static_cast<double>(i);
                const double posy = static_cast<double>(j);

                const double x = std::clamp((posx-dt*_grid.getU(i,j,b)), 0.0, static_cast<double>(F.x()+0.0));
                const double y = std::clamp((posy-dt*_grid.getV(i,j,b)), 0.0, static_cast<double>(F.y()+0.0));

                F(i,j) = interp(Fprev, x, y);
            }
        }
    }
}

void Fluids::preparePressureSolving(Laplacian& laplacian, Eigen::VectorXd& b)
{
    laplacian.A = Eigen::SparseMatrix<double>(_grid.activeCellsNb(), _grid.activeCellsNb());
    laplacian.A.reserve(Eigen::VectorXi::Constant(_grid.activeCellsNb(), 5));
    _grid._Adiag.reset();
    _grid._Ax.reset();
    _grid._Ay.reset();
    _grid._precon.reset();
    _grid._Az.reset();

    const double scale = 1;
    for (std::uint64_t n = 0; n < _grid._surface.maxIt(); ++n)
    {
        const std::uint64_t m = n % _grid._surface.maxIt();
        const std::uint16_t i = m % _grid._surface.x();
        const std::uint16_t j = m / _grid._surface.x();
        std::uint64_t id = _grid._pressureID(i,j);
        if (id > 0)
        {
            const std::uint64_t realID = id-1;
            std::uint64_t neibID = 0;
            if (i > 0 && (neibID = _grid._pressureID(i-1,j)) > 0)
            {
                laplacian.A.coeffRef(realID, neibID-1) = -scale;
                _grid._Adiag(i,j) += scale;
            }
            if (i+1 < Config::N && (neibID = _grid._pressureID(i+1,j)) > 0)
            {
                laplacian.A.coeffRef(realID, neibID-1) = -scale;
                _grid._Adiag(i,j) += scale;
                _grid._Ax(i,j) = -scale;
            }
            else if (i+1 < Config::N && _grid._surface.label(i+1,j) & EMPTY)
            {
                _grid._Adiag(i,j) += scale;
            }
            if (j > 0 && (neibID = _grid._pressureID(i,j-1)) > 0)
            {
                laplacian.A.coeffRef(realID, neibID-1) = -scale;
                _grid._Adiag(i,j) += scale;
            }
            if (j+1 < Config::N && (neibID = _grid._pressureID(i,j+1)) > 0)
            {
                laplacian.A.coeffRef(realID, neibID-1) = -scale;
                _grid._Adiag(i,j) += scale;
                _grid._Ay(i,j) = -scale;
            }
            else if (j+1 < Config::N && _grid._surface.label(i,j+1) & EMPTY)
            {
                _grid._Adiag(i,j) += scale;
            }

            std::uint8_t nonSolidNeib = 4;
            if (i == 0)
                nonSolidNeib--;
            if (j == 0)
                nonSolidNeib--;
            if (i == _grid._surface.x()-1)
                nonSolidNeib--;
            if (j == _grid._surface.y()-1)
                nonSolidNeib--;
            laplacian.A.coeffRef(realID, realID) = nonSolidNeib;
            b.coeffRef(realID) = div(i,j);
        }
    }
}

inline double Fluids::div(const std::uint16_t i, const std::uint16_t j) const
{
    const double h = 1.0/Config::N;
    double Udiv = _grid._U(i+1,j) - _grid._U(i,j);
    double Vdiv = _grid._V(i,j+1) - _grid._V(i,j);
    double rhs = - h * (Udiv + Vdiv);
    return rhs;
}

inline double Fluids::pressureAt(const std::uint16_t i, const std::uint16_t j, const Eigen::VectorXd x, const std::uint64_t l) const
{
    return _grid._pressure(i,j);
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
        for (std::uint16_t j = 0; j < std::max(_grid._U.y(), _grid._V.y()); ++j)
        {
            for (std::uint16_t i = 0; i < std::max(_grid._U.x(), _grid._V.x()); ++i)
            {
                if (j < _grid._surface.y() && i < _grid._surface.x() && (id = _grid._pressureID(i,j)) > 0)
                {
                    _grid._pressure(i,j) = x(id-1);
                }
                if (j < _grid._U.y() && i < _grid._U.x() && _grid._U.label(i,j) == LIQUID)
                {
                    _grid._U(i,j) -= Config::N*(pressureAt(i,j,x,0) - pressureAt(i-1,j,x,0));
                }
                if (j < _grid._V.y() && i < _grid._V.x() && _grid._V.label(i,j) == LIQUID)
                {
                    _grid._V(i,j) -= Config::N*(pressureAt(i,j,x,0) - pressureAt(i,j-1,x,0));
                }
            }
        }
    }
}

void Fluids::updateTexture()
{
    std::uint64_t it = 0;
    //double sum = 0;
    for (std::uint16_t i = 0; i < Config::N; ++i)
    {
        for (std::uint16_t j = 0; j < Config::N; ++j)
        {
            /*
            double value = _substance(i,j);
		    texture[it*3] = static_cast<std::uint8_t>(std::clamp(value, 0.0, 255.0));
		    texture[it*3+1] = static_cast<std::uint8_t>(std::clamp(value, 0.0, 255.0));
		    texture[it*3+2] = static_cast<std::uint8_t>(std::clamp(value, 0.0, 255.0));
            sum += value;
            */
            //const double implicit = _grid._surface(i,j);
            const double implicit = _grid._surface(i,j) <= 0 ? -255 : 255;
            double p = 0;
            double pneg = 0;
            if (_grid._pressure(i,j) < 0)
            {
                //pneg = _grid._pressure(i,j)*(-1);
            }
            else
            {
                //p = _grid._pressure(i,j);
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
                //_texture[it*3+0] = std::clamp(implicit*50, 0.0, 255.0);
                //_texture[it*3+1] = std::clamp(implicit*50, 0.0, 255.0);
                //_texture[it*3+2] = std::clamp(implicit*50, 0.0, 255.0);
            }
            /*
            if (_grid.activeCellsNb().find(_grid.hash(i,j)) != _grid.activeCellsNb().end())
            {
                _texture[it*3+0] = std::clamp(_texture[it*3+0] + 100.0, 0.0, 255.0);
            }
            */

            /*
            const double gO = _grid._surface.gradLength(i, j);
            _texture[it*3] = std::clamp(gO*96, 0.0, 255.0);
            _texture[it*3+1] = 0;
            _texture[it*3+2] = 0;
            */

            it++;
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

const std::uint16_t& Fluids::N() const
{
    return Config::N;
}

bool Fluids::isCellActive(std::uint16_t i, std::uint16_t j) const
{
    return _grid._pressureID(i,j) > 0;
}
