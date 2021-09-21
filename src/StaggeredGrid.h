#pragma once

#include <iomanip>
#include <algorithm>
#include <vector>

#include "./Eigen/Sparse"
#include "./utils.h"
#include "./config.h"

enum CellLabel
{
    EMPTY           = (1 << 0),
    LIQUID          = (1 << 1),
    SOLID           = (1 << 2),
    EXTRAPOLATED    = (1 << 3),
};

constexpr inline CellLabel operator|(CellLabel a, CellLabel b)
{
    return static_cast<CellLabel>(static_cast<int>(a) | static_cast<int>(b));
}

template<typename T, typename U>
class Field
{
 public:
    explicit Field(U Xsize, U Ysize, U Zsize)
        : _Xsize(Xsize)
        , _Ysize(Ysize)
        , _Zsize(Zsize)
    {
        if (Config::dim == 2)
        {
            _Zsize = 1;
        }
        _grid.resize(_Xsize*_Ysize*_Zsize);
        _label.resize(_Xsize*_Ysize*_Zsize);
        for (std::uint64_t it = 0; it < _Xsize*_Ysize*_Zsize; ++it)
        {
            _grid[it] = 0;
            _label[it] = EMPTY;
        }
        _maxIt = _Xsize*_Ysize*_Zsize;
    }

    T& operator()(const std::uint64_t idx)
    {
        return _grid[idx];
    }
    const T& operator()(const std::uint64_t idx) const
    {
        return _grid[idx];
    }
    T& operator()(const U i, const U j, const U k)
    {
        return _grid[idx(i, j, k)];
    }
    const T& operator()(const U i, const U j, const U k) const
    {
        return _grid[idx(i, j, k)];
    }
    const U& x() const
    {
        return _Xsize;
    }
    const U& y() const
    {
        return _Ysize;
    }
    const U& z() const
    {
        return _Zsize;
    }
    const std::uint64_t& maxIt() const
    {
        return _maxIt;
    }
    const std::vector<T>& data() const
    {
        return _grid;
    }
    const CellLabel& label(const U i, const U j, const U k) const
    {
        return _label[idx(i, j, k)];
    }
    CellLabel& label(const U i, const U j, const U k)
    {
        return _label[idx(i, j, k)];
    }
    CellLabel& label(const std::uint64_t idx)
    {
        return _label[idx];
    }
    bool checked(const U i, const U j, const U k) const
    {
        return _label[idx(i, j, k)] & LIQUID
            || _label[idx(i, j, k)] & EXTRAPOLATED;
    }
    void reset()
    {
        std::fill(_grid.begin(), _grid.end(), 0.0);
        for (std::uint64_t it = 0; it < _maxIt; ++it)
        {
            label(it) = EMPTY;
        }
    }
    void resetPos()
    {
        for (std::uint64_t it = 0; it < _maxIt; ++it)
        {
            label(it) = EMPTY;
        }
    }
    void setLabels(Field<T, U>& FU, Field<T, U>& FV, Field<T, U>& FW)
    {
        FU.resetPos();
        FV.resetPos();
        FW.resetPos();
        for (std::uint16_t k = 0; k < z(); ++k)
        {
            for (std::uint16_t j = 0; j < y(); ++j)
            {
                for (std::uint16_t i = 0; i < x(); ++i)
                {
                    if (operator()(i, j, k) < 0.0)
                    {
                        FU.label(i, j, k) = LIQUID;
                        FU.label(i+1, j, k) = LIQUID;
                        FV.label(i, j, k) = LIQUID;
                        FV.label(i, j+1, k) = LIQUID;
                        FW.label(i, j, k) = LIQUID;
                        FW.label(i, j, k+1) = LIQUID;
                    }
                }
            }
        }
        for (std::uint16_t k = 0;
                k < std::max({FU.z(), FV.z(), FW.z()});
                ++k)
        {
            for (std::uint16_t j = 0;
                    j < std::max({FU.y(), FV.y(), FW.y()});
                    ++j)
            {
                for (std::uint16_t i = 0;
                        i < std::max({FU.x(), FV.x(), FW.x()});
                        ++i)
                {
                    if (k < FU.z() && j < FU.y() && i < FU.x())
                    {
                        if (i == 0 || i == FU.x()-1)
                        {
                            FU.label(i, j, k) = SOLID;
                        }
                    }
                    if (k < FV.z() && j < FV.y() && i < FV.x())
                    {
                        if (j == 0 || j == FV.y()-1)
                        {
                            FV.label(i, j, k) = SOLID;
                        }
                    }
                    if (k < FW.z() && j < FW.y() && i < FW.x())
                    {
                        if (k == 0 || k == FW.z()-1)
                        {
                            FW.label(i, j, k) = SOLID;
                        }
                    }
                }
            }
        }
    }
    void setFromVec(const Eigen::VectorXd& v)
    {
        for (std::uint64_t it = 0; it < _maxIt; ++it)
        {
            _grid[it] = v.coeff(it);
        }
    }
    Eigen::VectorXd vec() const
    {
        Eigen::VectorXd v(_maxIt);
        for (std::uint64_t it = 0; it < _maxIt; ++it)
        {
            v.coeffRef(it) = _grid[it];
        }
        return v;
    }
    inline double gradLength(const U i, const U j, const U k) const
    {
        double gradI = 0.0;
        double gradJ = 0.0;
        double gradK = 0.0;
        if (i == 0)
        {
            gradI = operator()(0, j, k) - operator()(1, j, k);
        }
        else if (i == _Xsize-1)
        {
            gradI = operator()(_Xsize-2, j, k) - operator()(_Xsize-1, j, k);
        }
        else
        {
            if (std::abs(operator()(i+1, j, k))
                    < std::abs(operator()(i-1, j, k)))
            {
                gradI = operator()(i, j, k) - operator()(i+1, j, k);
            }
            else
            {
                gradI = operator()(i-1, j, k) - operator()(i, j, k);
            }
        }
        if (j == 0)
        {
            gradJ = operator()(i, 0, k) - operator()(i, 1, k);
        }
        else if (j == _Ysize-1)
        {
            gradJ = operator()(i, _Ysize-2, k) - operator()(i, _Ysize-1, k);
        }
        else
        {
            if (std::abs(operator()(i, j+1, k))
                    < std::abs(operator()(i, j-1, k)))
            {
                gradJ = operator()(i, j, k) - operator()(i, j+1, k);
            }
            else
            {
                gradJ = operator()(i, j-1, k) - operator()(i, j, k);
            }
        }
        if (k == 0)
        {
            gradK = operator()(i, j, 0) - operator()(i, j, 1);
        }
        else if (k == _Zsize-1)
        {
            gradK = operator()(i, j, _Zsize-2) - operator()(i, j, _Zsize-1);
        }
        else
        {
            if (std::abs(operator()(i, j, k+1))
                    < std::abs(operator()(i, j, k-1)))
            {
                gradK = operator()(i, j, k) - operator()(i, j, k+1);
            }
            else
            {
                gradK = operator()(i, j, k-1) - operator()(i, j, k);
            }
        }
        return std::sqrt(
                  std::pow(gradI, 2)
                + std::pow(gradJ, 2)
                + std::pow(gradK, 2)
        );
    }

    friend std::ostream& operator<<(std::ostream& os, const Field& obj)
    {
        os << std::fixed << std::setprecision(2);
        for (std::int64_t j = obj._Ysize-1; j >= 0; --j)
        {
            for (std::int64_t k = 0; k < obj._Zsize; ++k)
            {
                for (std::int64_t i = 0; i < obj._Xsize; ++i)
                {
                    if (obj(i, j, k) >= 0.0)
                    {
                        os << " ";
                    }
                    os << obj(i, j, k) << " ";
                }
                os << std::endl;
            }
            os << std::endl << std::endl;
        }
        return os;
    }

    inline std::uint64_t idx(const U i, const U j, const U k) const
    {
        if (i > _Xsize)
        {
            std::cout << i << std::endl;
            ERROR("i > _Xsize");
        }
        if (j > _Ysize)
        {
            std::cout << j << std::endl;
            ERROR("j > _Ysize");
        }
        if (k > _Zsize)
        {
            std::cout << k << std::endl;
            ERROR("k > _Zsize");
        }
        return i + j * _Xsize + k * _Xsize * _Ysize;
    }

 private:
    std::uint64_t _maxIt;
    std::vector<T> _grid;
    std::vector<CellLabel> _label;
    U _Xsize;
    U _Ysize;
    U _Zsize;
};

template<typename T, typename R>
class StaggeredGrid
{
 public:
    explicit StaggeredGrid(R N) : _N(N) {}
    inline std::uint64_t hash(
            const std::uint16_t i,
            const std::uint16_t j,
            const std::uint16_t k
        ) const
    {
        return (static_cast<std::uint64_t>(i) << 32) | (j << 16) | k;
    }
    inline const T getU(
            const R i,
            const R j,
            const R k,
            const std::uint8_t b
        ) const
    {
        switch (b)
        {
            case 0:
                return 0.5*(_U(i, j, k)+_U(i+1, j, k));
                break;
            case 1:
                return _U(i, j, k);
                break;
            case 2:
                if (j == 0)
                {
                    return 0.5*(_U(i, j, k)+_U(i+1, j, k));
                }
                else if (j == _N-1)
                {
                    return 0.5*(_U(i, j-1, k)+_U(i+1, j-1, k));
                }
                else
                {
                    return 0.25*(_U(i, j, k)+_U(i+1, j, k)+_U(i, j-1, k)+_U(i+1, j-1, k));
                }
                break;
            case 3:
                if (k == 0)
                {
                    return 0.5*(_U(i, j, k)+_U(i+1, j, k));
                }
                else if (k == _N-1)
                {
                    return 0.5*(_U(i, j, k-1)+_U(i+1, j, k-1));
                }
                else
                {
                    return 0.25*(_U(i, j, k)+_U(i+1, j, k)+_U(i, j, k-1)+_U(i+1, j, k-1));
                }
                break;
        }
        ERROR("b should be either 0, 1, 2 or 3");
    }
    inline const T getV(
            const R i,
            const R j,
            const R k,
            const std::uint8_t b
        ) const
    {
        switch (b)
        {
            case 0:
                return 0.5*(_V(i, j, k)+_V(i, j+1, k));
                break;
            case 1:
                if (i == 0)
                {
                    return 0.5*(_V(i, j, k)+_V(i, j+1, k));
                }
                else if (i == _N-1)
                {
                    return 0.5*(_V(i-1, j, k)+_V(i-1, j+1, k));
                }
                else
                {
                    return 0.25*(_V(i, j, k)+_V(i, j+1, k)+_V(i-1, j, k)+_V(i-1, j+1, k));
                }
                break;
            case 2:
                return _V(i, j, k);
                break;
            case 3:
                if (k == 0)
                {
                    return 0.5*(_V(i, j, k)+_V(i, j+1, k));
                }
                else if (k == _N-1)
                {
                    return 0.5*(_V(i, j, k-1)+_V(i, j+1, k-1));
                }
                else
                {
                    return 0.25*(_V(i, j, k)+_V(i, j+1, k)+_V(i, j, k-1)+_V(i, j+1, k-1));
                }
                break;
        }
        ERROR("b should be either 0, 1, 2 or 3");
    }
    inline const T getW(
            const R i,
            const R j,
            const R k,
            const std::uint8_t b
        ) const
    {
        switch (b)
        {
            case 0:
                return 0.5*(_W(i, j, k)+_W(i, j, k+1));
                break;
            case 1:
                if (i == 0)
                {
                    return 0.5*(_W(i, j, k)+_W(i, j, k+1));
                }
                else if (i == _N-1)
                {
                    return 0.5*(_W(i-1, j, k)+_W(i-1, j, k+1));
                }
                else
                {
                    return 0.25*(_W(i, j, k)+_W(i, j, k+1)+_W(i-1, j, k)+_W(i-1, j, k+1));
                }
                break;
            case 2:
                if (j == 0)
                {
                    return 0.5*(_W(i, j, k)+_W(i, j, k+1));
                }
                else if (j == _N-1)
                {
                    return 0.5*(_W(i, j-1, k)+_W(i, j-1, k+1));
                }
                else
                {
                    return 0.25*(_W(i, j, k)+_W(i, j, k+1)+_W(i, j-1, k)+_W(i, j-1, k+1));
                }
                break;
            case 3:
                return _W(i, j, k);
                break;
        }
        ERROR("b should be either 0, 1, 2 or 3");
    }
    void tagActiveCells()
    {
        _activeCells = 0;
        _pressureID.reset();
        std::uint64_t id = 0;
        for (std::uint16_t k = 0; k < _surface.z(); ++k)
        {
            for (std::uint16_t j = 0; j < _surface.y(); ++j)
            {
                for (std::uint16_t i = 0; i < _surface.x(); ++i)
                {
                    if (_surface(i, j, k) < 0.0)
                    {
                        _surface.label(i, j, k) = LIQUID;
                        _pressureID(i, j, k) = ++id;
                        _activeCells++;
                    }
                    else
                    {
                        _surface.label(i, j, k) = EMPTY;
                    }
                }
            }
        }
    }
    inline std::uint64_t activeCellsNb()
    {
        return _activeCells;
    }

    R _N;
    Field<T, R> _substance {_N, _N, _N};
    Field<T, R> _surface {_N, _N, _N};
    Field<T, R> _U {static_cast<std::uint16_t>(_N+1), _N, _N};
    Field<T, R> _V {_N, static_cast<std::uint16_t>(_N+1), _N};
    Field<T, R> _W {_N, _N, static_cast<std::uint16_t>(_N+1)};
    Field<T, R> _UPrev {static_cast<std::uint16_t>(_N+1), _N, _N};
    Field<T, R> _VPrev {_N, static_cast<std::uint16_t>(_N+1), _N};
    Field<T, R> _WPrev {_N, _N, static_cast<std::uint16_t>(_N+1)};
    Field<T, R> _pressure {_N, _N, _N};
    Field<T, R> _Adiag {_N, _N, _N};
    Field<T, R> _Ax {_N, _N, _N};
    Field<T, R> _Ay {_N, _N, _N};
    Field<T, R> _Az {_N, _N, _N};
    Field<T, R> _precon {_N, _N, _N};
    Field<T, R> _q {_N, _N, _N};
    Field<T, R> _z {_N, _N, _N};
    Field<std::uint64_t, R> _pressureID {_N, _N, _N};
    Field<T, R> _substancePrev {_N, _N, _N};
    Field<T, R> _surfacePrev {_N, _N, _N};

 private:
    std::uint64_t _activeCells {0};
};
