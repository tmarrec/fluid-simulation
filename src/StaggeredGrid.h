#pragma once

#include "utils.h"

#include <Eigen/Sparse>
#include <iomanip>

template<typename T, typename U>
class Field
{
public:
    Field(U Xsize, U Ysize) : _Xsize(Xsize), _Ysize(Ysize)
    {
        _grid.resize(Xsize*Ysize);
        _gridSet.resize(Xsize*Ysize);
        for (std::uint64_t it = 0; it < Xsize*Ysize; ++it)
        {
            _grid[it] = 0;
            _gridSet[it] = false;
        }
        _maxIt = Xsize*Ysize;
    }
            T& operator()(const U idx)                                  { return _grid[idx]; }
    const   T& operator()(const U idx)                          const   { return _grid[idx]; }
            T& operator()(const U i, const U j)                         { return _grid[idx(i,j)]; }
    const   T& operator()(const U i, const U j)                 const   { return _grid[idx(i,j)]; }
    const   U& x()                                              const   { return _Xsize; }
    const   U& y()                                              const   { return _Ysize; }
    const   std::uint64_t& maxIt()                              const   { return _maxIt; }
    const   std::vector<T>& data()                              const   { return _grid; }
            bool isSet(const U i, const U j)                    const   { return _gridSet[idx(i,j)]; }
            void set(const U i, const U j, const bool value)            { _gridSet[idx(i,j)] = value; }
            void set(const U idx, const bool value)                     { _gridSet[idx] = value; }
            void reset()
            {
                std::fill(_grid.begin(), _grid.end(), 0.0);
                for (std::uint64_t it = 0; it < _maxIt; ++it)
                {
                    set(it, false);
                }
            };
            void resetBool()
            {
                for (std::uint64_t it = 0; it < _maxIt; ++it)
                {
                    set(it, false);
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
    std::vector<bool> _gridSet;
    U _Xsize;
    U _Ysize;
};

template<typename T, typename U>
class StaggeredGrid
{
public:
    StaggeredGrid(U N) : _N(N) {}
    inline const T getU(const U i, const U j, const std::uint8_t b) const
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

    inline const T getV(const U i, const U j, const std::uint8_t b) const
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

    U _N;
    Field<T, U> _substance {_N, _N};
    Field<T, U> _surface {_N, _N};
    Field<T, U> _U {static_cast<std::uint16_t>(_N+1), _N};
    Field<T, U> _V {_N, static_cast<std::uint16_t>(_N+1)};

    Field<T, U> _substancePrev {_N, _N};
    Field<T, U> _surfacePrev {_N, _N};
    Field<T, U> _UPrev {static_cast<std::uint16_t>(_N+1), _N};
    Field<T, U> _VPrev {_N, static_cast<std::uint16_t>(_N+1)};
private:
};
