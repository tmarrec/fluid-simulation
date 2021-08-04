#pragma once

#include "utils.h"

#include <Eigen/Sparse>
#include <iomanip>

enum CellPosition
{
    UNKNOWN,
    OUTSIDE,
    OUTSIDE_EXTRAPOLATED,
    INSIDE,
};

struct Cell
{
    std::uint16_t i;
    std::uint16_t j;
    std::uint64_t label;
};

template<typename T, typename U>
class Field
{
public:
    explicit Field(U Xsize, U Ysize) : _Xsize(Xsize), _Ysize(Ysize)
    {
        _grid.resize(Xsize*Ysize);
        _pos.resize(Xsize*Ysize);
        for (std::uint64_t it = 0; it < Xsize*Ysize; ++it)
        {
            _grid[it] = 0;
            _pos[it] = UNKNOWN;
        }
        _maxIt = Xsize*Ysize;
    }
            T& operator()(const std::uint64_t idx)                      { return _grid[idx]; }
    const   T& operator()(const std::uint64_t idx)              const   { return _grid[idx]; }
            T& operator()(const U i, const U j)                         { return _grid[idx(i,j)]; }
    const   T& operator()(const U i, const U j)                 const   { return _grid[idx(i,j)]; }
    const   U& x()                                              const   { return _Xsize; }
    const   U& y()                                              const   { return _Ysize; }
    const   std::uint64_t& maxIt()                              const   { return _maxIt; }
    const   std::vector<T>& data()                              const   { return _grid; }
    const   CellPosition& pos(const U i, const U j)             const   { return _pos[idx(i,j)]; }
            CellPosition& pos(const U i, const U j)                     { return _pos[idx(i,j)]; }
            CellPosition& pos(const std::uint64_t idx)                  { return _pos[idx]; }
            bool checked(const U i, const U j)                  const   { return _pos[idx(i,j)] == INSIDE || _pos[idx(i,j)] == OUTSIDE_EXTRAPOLATED; }
            void reset()
            {
                std::fill(_grid.begin(), _grid.end(), 0.0);
                for (std::uint64_t it = 0; it < _maxIt; ++it)
                {
                    pos(it) = UNKNOWN;
                }
            };
            void resetPos()
            {
                for (std::uint64_t it = 0; it < _maxIt; ++it)
                {
                    pos(it) = UNKNOWN;
                }
            };
            void setFromVec(const Eigen::VectorXd& v)
            {
                for (std::uint64_t it = 0; it < _maxIt; ++it)
                {
                    _grid[it] = v.coeff(it);
                }
            }
            Eigen::VectorXd vec()                               const
            { 
                Eigen::VectorXd v(_maxIt);
                for (std::uint64_t it = 0; it < _maxIt; ++it)
                {
                    v.coeffRef(it) = _grid[it];
                }
                return v;
            }
            inline double gradLength(const U i, const U j)      const
            {
                double gradI = 0;
                double gradJ = 0;
                if (i == 0)
                {
                    gradI = operator()(0,j) - operator()(1,j);
                }
                else if (i == _Xsize-1)
                {
                    gradI = operator()(_Xsize-2,j) - operator()(_Xsize-1,j);
                }
                else
                {
                    if (std::abs(operator()(i+1,j)) < std::abs(operator()(i-1,j)))
                    {
                        gradI = operator()(i,j) - operator()(i+1,j);
                    }
                    else
                    {
                        gradI = operator()(i-1,j) - operator()(i,j);
                    }
                }
                if (j == 0)
                {
                    gradJ = operator()(i,0) - operator()(i,1);
                }
                else if (j == _Ysize-1)
                {
                    gradJ = operator()(i,_Ysize-2) - operator()(i,_Ysize-1);
                }
                else
                {
                    if (std::abs(operator()(i,j+1)) < std::abs(operator()(i,j-1)))
                    {
                        gradJ = operator()(i,j) - operator()(i,j+1);
                    }
                    else
                    {
                        gradJ = operator()(i,j-1) - operator()(i,j);
                    }
                }
                return std::sqrt(std::pow(gradI, 2)+std::pow(gradJ, 2));
            }

    friend std::ostream& operator<<(std::ostream& os, const Field& obj)
    {
        os << std::fixed << std::setprecision(2);
        for (U j = 0; j < obj._Ysize; ++j)
        {
            for (U i = 0; i < obj._Xsize; ++i)
            {
                if (obj(i,j) >= 0)
                {
                    os << " ";
                }
                os << obj(i,j) << " ";
            }
            os << std::endl;
        }
        return os;
    }

    inline std::uint64_t idx(const U i, const U j) const
    { 
        return i + j * _Xsize;
    };

private:
    std::uint64_t _maxIt;
    std::vector<T> _grid;
    std::vector<CellPosition> _pos;
    U _Xsize;
    U _Ysize;
};

template<typename T, typename R>
class StaggeredGrid
{
public:
    explicit StaggeredGrid(R N) : _N(N) {}
    inline std::uint64_t hash(const std::uint16_t i, const std::uint16_t j) const
    {
        return (i << 16) | j;
    }
    inline const T getU(const R i, const R j, const std::uint8_t b) const
    {
        switch (b)
        {
            case 0:
                return 0.5*(_U(i,j)+_U(i+1,j));
                break;
            case 1:
                return _U(i,j);
                break;
            case 2:
                if (j == 0)
                {
                    return 0.5*(_U(i,j)+_U(i+1,j));
                }
                else if (j == _N-1)
                {
                    return 0.5*(_U(i,j-1)+_U(i+1,j-1));
                }
                else
                {
                    return 0.25*(_U(i,j)+_U(i+1,j)+_U(i,j-1)+_U(i+1,j-1));
                }
                break;

        }
        ERROR("b should be either 0, 1 or 2");
    }
    inline const T getV(const R i, const R j, const std::uint8_t b) const
    {
        switch (b)
        {
            case 0:
                return 0.5*(_V(i,j)+_V(i,j+1));
                break;
            case 1:
                if (i == 0)
                {
                    return 0.5*(_V(i,j)+_V(i,j+1));
                }
                else if (i == _N-1)
                {
                    return 0.5*(_V(i-1,j)+_V(i-1,j+1));
                }
                else
                {
                    return 0.25*(_V(i,j)+_V(i,j+1)+_V(i-1,j)+_V(i-1,j+1));
                }
                break;
            case 2:
                return _V(i,j);
                break;

        }
        ERROR("b should be either 0, 1 or 2");
    }
    void tagActiveCells()
    {
        _activeCells.clear();
        std::uint64_t label = 1;
        for (std::uint16_t j = 0; j < _surface.y(); ++j)
        {
            for (std::uint16_t i = 0; i < _surface.x(); ++i)
            {
                if (_surface(i,j) < 0.0)
                {
                    _activeCells.insert({hash(i,j), {i,j,label}});
                    label++;
                }
            }
        }
    }

    R _N;
    Field<T, R> _substance {_N, _N};
    Field<T, R> _surface {_N, _N};
    Field<T, R> _U {static_cast<std::uint16_t>(_N+1), _N};
    Field<T, R> _V {_N, static_cast<std::uint16_t>(_N+1)};
    Field<T, R> _pressure {_N, _N};

    Field<T, R> _substancePrev {_N, _N};
    Field<T, R> _surfacePrev {_N, _N};
    Field<T, R> _UPrev {static_cast<std::uint16_t>(_N+1), _N};
    Field<T, R> _VPrev {_N, static_cast<std::uint16_t>(_N+1)};

    std::unordered_map<std::uint64_t, Cell> _activeCells;
private:
};
