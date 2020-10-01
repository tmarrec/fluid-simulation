#include "BSpline.h"

#include <iostream>

BSpline::BSpline(
		unsigned short order,
		std::vector<float> knots,
		std::vector<glm::vec3> controls
	)
	: _order(order)
	, _knots(knots)
	, _controls(controls)
{
	_range = std::tuple<float, float>(knots[order-1], knots[controls.size()]);
}

std::tuple<float, float> BSpline::range() const {
	return _range;
}

glm::vec3 BSpline::eval(float u) {
	unsigned short dec = 0;
	unsigned short i = _order;
	unsigned short m = _order-1;

	while (u > _knots[i]) {
		i++;
		dec++;
	}

	std::vector<glm::vec3> P_temps;
	unsigned short k = _order;
	
	for(unsigned a = dec; a < dec+_order; ++a)
		P_temps.push_back(_controls[a]);
	
	for (unsigned short l = 0; l < m; ++l) {
		for (unsigned short j = 0; j < k-1; ++j) {
			P_temps[j] = ((_knots[dec+k+j]-u)/(_knots[dec+k+j]-_knots[dec+1+j]))*P_temps[j]
						 +
						 ((u-_knots[dec+1+j])/(_knots[dec+k+j]-_knots[dec+1+j]))*P_temps[j+1];
		}
		dec++;
		k--;
	}

	return P_temps[0];
}

BSpline::~BSpline() {
	
}
