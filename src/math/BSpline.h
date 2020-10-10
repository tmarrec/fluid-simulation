#pragma once

#include <cstdint>
#include <vector>
#include <tuple>

#include "../utils.h"

struct Range
{
	double start;
	double end;
};

class BSpline
{
public:
	BSpline(std::uint8_t __order, std::vector<float> __knots, std::vector<glm::vec3> __controls)
	: _order { __order }
	, _knots { __knots }
	, _controls { __controls }
	, _range { __knots[__order-1], __knots[__controls.size()] }
	{}


	glm::vec3 eval(float u) const
	{
		std::uint8_t dec = 0;
		std::uint8_t i = _order;
		std::uint8_t m = _order-1;
		std::uint8_t k = _order;
		std::vector<glm::vec3> P_temps;


		while (u > _knots[i]) {
			i++;
			dec++;
		}
	
		for(unsigned a = dec; a < dec+_order; ++a)
			P_temps.push_back(_controls[a]);
	
		for (std::uint8_t l = 0; l < m; ++l) {
			for (std::uint8_t j = 0; j < k-1; ++j) {
				P_temps[j] = ((_knots[dec+k+j]-u)/(_knots[dec+k+j]-_knots[dec+1+j]))*P_temps[j]
							 +
							 ((u-_knots[dec+1+j])/(_knots[dec+k+j]-_knots[dec+1+j]))*P_temps[j+1];
			}
			dec++;
			k--;
		}

		return P_temps.front();
	}

	Range range() const { return _range; }

	~BSpline() {};

private:
	std::uint8_t _order;
	std::vector<float> _knots;
	std::vector<glm::vec3> _controls;
	Range _range;
};

