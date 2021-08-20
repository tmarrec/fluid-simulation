#include "Fluids.h"
#include "ConjugateGradient.h"
#include <iomanip>

Fluids::Fluids()
{
    //initLaplacians(_N, _dt, _viscosity, _diffusion, _laplacianViscosityX, _laplacianViscosityY, _laplacianProject, _laplacianDiffuse);
}

void Fluids::update([[maybe_unused]] std::uint64_t iteration)
{
    _iteration = iteration;
    
    /* Testings */
    std::uint16_t N = _N/2;

    // LEVEL-SET Init
    double r = 25;
    if (iteration == 0)
    {
        glm::vec2 pt;
        pt.x = N/*-(N/2.1)*/;
        pt.y = N;
        particles.emplace_back(pt);
        pt.x = N+(N/2.1);
        pt.y = N;
        //particles.emplace_back(pt);

        for (std::uint16_t j = 0; j < _N; ++j)
        {
            for (std::uint16_t i = 0; i < _N; ++i)
            {
                _grid._surface(i,j) = 10;
            }
        }

        /*
        for (std::uint16_t j = _N; j > _N-16; --j)
        {
            for (std::uint16_t i = 0; i < _N; ++i)
            {
                _grid._surface(i,j) = -10;
            }
        }
        */
    }

    if (_iteration == 0)
    {
        for (std::uint16_t j = 0; j < _N; ++j)
        {
            for (std::uint16_t i = 0; i < _N; ++i)
            {
                if (particles.size() > 0)
                {
                    double dist = std::sqrt(std::pow(i-particles.front().x,2)+std::pow(j-particles.front().y,2))-r;
                    for (const auto& p : particles)
                    {
                        dist = std::min(dist, std::sqrt(std::pow(i-p.x,2)+std::pow(j-p.y,2))-r);
                    }
                    if (dist <= 0)
                    _grid._surface(i,j) = dist;
                }
            }
        }
    }

    step();

    updateTexture();

    //std::cin.get();
    usleep(1000);
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
                _grid._U.pos(i,j) = LIQUID;
                _grid._U.pos(i+1,j) = LIQUID;
                _grid._V.pos(i,j) = LIQUID;
                _grid._V.pos(i,j+1) = LIQUID;
            }
        }
    }
    for (std::uint16_t j = 0; j < _grid._surface.y()+1; ++j)
    {
        for (std::uint16_t i = 0; i < _grid._surface.x()+1; ++i)
        {
            if (j < _grid._U.y() && (i == 0 || i == _grid._surface.x()))
            {
                _grid._U.pos(i,j) = SOLID;
            }
            if (i < _grid._V.x() && (j == 0 || j == _grid._surface.y()))
            {
                _grid._V.pos(i,j) = SOLID;
            }
        }
    }

    // Extrapolate the velocity field
    extrapolate(_grid._U);
    extrapolate(_grid._V);

    // Advect level-set everywhere using the fully extrapolated velocity
    advect(_grid._surface, _grid._surfacePrev, 0);
    redistancing(3, _grid._surface);

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
    const double G = 0.1;
    for (std::uint16_t j = 0; j < _grid._V.y(); ++j)
    {
        for (std::uint16_t i = 0; i < _grid._V.x(); ++i)
        {
            if (!(_grid._V.pos(i,j) & SOLID))
            {
                _grid._V(i,j) += G;
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

void Fluids::extrapolate(Field<double,std::uint16_t>& F, std::uint16_t nbIte)
{
    std::uint16_t it = 0;
    Field<double,std::uint16_t> temp {F.x(), F.y()};
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

            if (F.pos(i,j) & LIQUID)
            {
                temp(i,j) = F(i,j);
                temp.pos(i,j) = LIQUID;
            }
            else if (F.pos(i,j) & SOLID)
            {
                temp(i,j) = F(i,j);
                temp.pos(i,j) = SOLID;
            }
            else if (F.pos(i,j) & EXTRAPOLATED)
            {
                temp(i,j) = F(i,j);
                temp.pos(i,j) = temp.pos(i,j) | EXTRAPOLATED;
            }
            else
            {
                std::uint8_t nbNeighbors = 0;
                double value = 0.0;
                if (i < _N-1 && F.checked(i+1,j))
                {
                    nbNeighbors++;
                    value += F(i+1,j);
                }
                if (i > 0 && F.checked(i-1,j))
                {
                    nbNeighbors++;
                    value += F(i-1,j);
                }
                if (j < _N-1 && F.checked(i,j+1))
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
                    temp(i,j) = value/nbNeighbors;
                    temp.pos(i,j) = temp.pos(i,j) | EXTRAPOLATED;
                }
            }
        }
        F = temp;
        if (nbIte > 0)
        {
            if (++it == nbIte)
            {
                return;
            }
        }
    } while (nbNeg > 0);
}

void Fluids::redistancing(const std::uint64_t nbIte, Field<double, std::uint16_t>& field)
{
    const double dx = 1.0/_N;
    auto F = field;
    auto QNew = field;
    for (std::uint16_t j = 0; j < field.y(); ++j)
    {
        for (std::uint16_t i = 0; i < field.x(); ++i)
        {
            if (    (i+1 < field.x() && !(field(i,j) >= 0) ^ (field(i+1,j) < 0)) ||
                    (i-1 >= 0 && !(field(i,j) >= 0) ^ (field(i-1,j) < 0)) ||
                    (j+1 < field.y() && !(field(i,j) >= 0) ^ (field(i,j+1) < 0)) ||
                    (j-1 >= 0 && !(field(i,j) >= 0) ^ (field(i,j-1) < 0))   )
            {
                F(i,j) = 1;
                F.pos(i,j) = LIQUID;
            }
            else
            {
                F(i,j) = 0;
                F.pos(i,j) = AIR;
            }
        }
    }
    extrapolate(F, 2);

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
    Eigen::VectorXd x(A.diag.size());
    ConjugateGradient(A, x, Fprev.vec(), _solverType);
    F.setFromVec(x);
}


void Fluids::advect(Field<double,std::uint16_t>& F, Field<double,std::uint16_t>& Fprev, const std::uint8_t b) const
{
    Fprev = F;
	const double dt = _dt * _N;
    for (std::uint16_t j = 0; j < F.y(); ++j)
    {
        for (std::uint16_t i = 0; i < F.x(); ++i)
        {
            if (!(F.pos(i,j) & SOLID))
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

void Fluids::pressureMatrix(Laplacian& laplacian, Eigen::VectorXd& b) const
{
    std::vector<Eigen::Triplet<double>> tripletListA;
    std::for_each(_grid._activeCells.begin(), _grid._activeCells.end(),
    [&](const auto& it)
    {
        const std::uint16_t i = it.second.i; 
        const std::uint16_t j = it.second.j; 
        const std::uint64_t label = it.second.label; 
        std::unordered_map<std::uint64_t,Cell>::const_iterator nCell;
        if (((nCell = _grid._activeCells.find(_grid.hash(i+1,j)))) != _grid._activeCells.end() && i+1 < _N)
        {
            tripletListA.emplace_back(Eigen::Triplet<double>(label-1, nCell->second.label-1, -1));
        }
        if (((nCell = _grid._activeCells.find(_grid.hash(i-1,j)))) != _grid._activeCells.end() && i != 0)
        {
            tripletListA.emplace_back(Eigen::Triplet<double>(label-1, nCell->second.label-1, -1));
        }
        if (((nCell = _grid._activeCells.find(_grid.hash(i,j+1)))) != _grid._activeCells.end() && j+1 < _N)
        {
            tripletListA.emplace_back(Eigen::Triplet<double>(label-1, nCell->second.label-1, -1));
        }
        if (((nCell = _grid._activeCells.find(_grid.hash(i,j-1)))) != _grid._activeCells.end() && j != 0)
        {
            tripletListA.emplace_back(Eigen::Triplet<double>(label-1, nCell->second.label-1, -1));
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
        tripletListA.emplace_back(Eigen::Triplet<double>(label-1, label-1, nonSolidNeib));

        b.coeffRef(label-1) = div(i,j);
    }); 

    Eigen::SparseMatrix<double> A = Eigen::SparseMatrix<double>(_grid._activeCells.size(), _grid._activeCells.size());

    A.setFromTriplets(tripletListA.begin(), tripletListA.end());
    laplacian.A = A;
    setAMatrices(laplacian);
    setPrecon(laplacian);
}

inline double Fluids::div(const std::uint16_t i, const std::uint16_t j) const
{
    const double h = 1.0/_N;
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
    if (_grid._activeCells.size() > 0)
    {
        Eigen::VectorXd x(_grid._activeCells.size());
        Eigen::VectorXd b(_grid._activeCells.size());

        Laplacian A {};
        pressureMatrix(A, b);
        ConjugateGradient(A, x, b, _solverType);

        _grid._pressure.reset();
        std::for_each(_grid._activeCells.begin(), _grid._activeCells.end(),
        [&](const auto& elem)
        {
            const std::uint16_t i = elem.second.i; 
            const std::uint16_t j = elem.second.j; 
            const std::uint64_t label = elem.second.label; 

            _grid._pressure(i,j) = x.coeff(label-1);
        });

        // Should optimize
        for (std::uint16_t j = 0; j < _grid._U.y(); ++j)
        {
            for (std::uint16_t i = 0; i < _grid._U.x(); ++i)
            {
                if (_grid._U.pos(i,j) == LIQUID)
                {
                    _grid._U(i,j) -= _N*(pressureAt(i,j,x,0) - pressureAt(i-1,j,x,0));
                }
            }
        }
        for (std::uint16_t j = 0; j < _grid._V.y(); ++j)
        {
            for (std::uint16_t i = 0; i < _grid._V.x(); ++i)
            {
                if (_grid._V.pos(i,j) == LIQUID)
                {
                    _grid._V(i,j) -= _N*(pressureAt(i,j,x,0) - pressureAt(i,j-1,x,0));
                }
            }
        }
    }
}

void Fluids::updateTexture()
{
    std::uint64_t it = 0;
    //double sum = 0;
    for (std::uint16_t i = 0; i < _N; ++i)
    {
        for (std::uint16_t j = 0; j < _N; ++j)
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
                pneg = _grid._pressure(i,j)*(-1);
            }
            else
            {
                p = _grid._pressure(i,j);
            }
            if (implicit <= 0)
            {
                _texture[it*3+0] = std::clamp(p*1500, 0.0, 255.0);
                _texture[it*3+1] = std::clamp(pneg*1500, 0.0, 255.0);
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
            if (_grid._activeCells.find(_grid.hash(i,j)) != _grid._activeCells.end())
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
    return _N;
}

bool Fluids::isCellActive(std::uint16_t i, std::uint16_t j) const
{
    return _grid._activeCells.find(_grid.hash(i,j)) != _grid._activeCells.end();
}

