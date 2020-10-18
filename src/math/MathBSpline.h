#pragma once

#include <cstdint>
#include <vector>
#include <tuple>

#include "../utils.h"

class MathBSpline
{
public:
	MathBSpline(std::uint8_t __order, std::vector<glm::vec3> __controls)
	: _order { __order }
	, _knots { _uniform_vector(__controls.size()+__order+1) }
	, _controls { __controls }
	, _range { _knots[__order-1], _knots[__controls.size()] }
	{}

	MathBSpline(std::uint8_t __order, std::vector<glm::vec3> __controls, std::vector<float> __knots)
	: _order { __order }
	, _knots { __knots }
	, _controls { __controls }
	, _range { _knots[__order-1], _knots[__controls.size()] }
	{}

	glm::vec3 eval(float __u) const
	{
		ASSERT(__u >= _range.start && __u <= _range.end, "__u is not in the range of the BSpline");
		std::uint8_t dec = 0;
		std::uint8_t i = _order;
		std::uint8_t m = _order-1;
		std::uint8_t k = _order;
		std::vector<glm::vec3> P_temps;
		while (__u > _knots[i])
		{
			i++;
			dec++;
		}
		for(unsigned a = dec; a < dec+_order; ++a)
		{
			P_temps.emplace_back(_controls[a]);
		}
		for (std::uint8_t l = 0; l < m; ++l)
		{
			for (std::uint8_t j = 0; j < k-1; ++j)
			{
				P_temps[j] = ((_knots[dec+k+j]-__u)/(_knots[dec+k+j]-_knots[dec+1+j]))*P_temps[j]+((__u-_knots[dec+1+j])/(_knots[dec+k+j]-_knots[dec+1+j]))*P_temps[j+1];
			}
			dec++;
			k--;
		}
		return P_temps.front();
	}

	Range range() const { return _range; }
	~MathBSpline() {};

private:
	std::vector<float> _uniform_vector(std::uint8_t __size) const
	{
		std::vector<float> v;
		for (std::uint8_t i = 0; i < __size; ++i)
		{
			v.emplace_back(i);
		}
		return v;
	}

	std::uint8_t _order;
	std::vector<float> _knots;
	std::vector<glm::vec3> _controls;
	Range _range;
};

