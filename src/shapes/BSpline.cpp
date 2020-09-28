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

	std::cout << "u: " << u << std::endl;
	std::cout << "_knots[i]: " << _knots[i] << std::endl;
	while (u > _knots[i]) {
		i++;
		dec++;
	}

	std::vector<glm::vec3> P_temps;
	unsigned short k = _order;
	
	for(unsigned a = dec; a < dec+_order; ++a)
		P_temps.push_back(_controls[a]);
	
	std::cout << "CONTROLS : ";
	for (auto p : _controls)
		std::cout << "(" << p.x << " " << p.y << " " << p.z << ") ";
	std::cout << std::endl;

	std::cout << "P_temps : ";
	for (auto p : P_temps)
		std::cout << "(" << p.x << " " << p.y << " " << p.z << ") ";
	std::cout << std::endl;

	for (unsigned short l = 0; l < m; ++l) {
		for (unsigned short j = 0; j < k-1; ++j) {
			P_temps[j] = ((_knots[dec+k+j]-u)/(_knots[dec+k+j]-_knots[dec+1+j]))*P_temps[j]
						 +
						 ((u-_knots[dec+1+j])/(_knots[dec+k+j]-_knots[dec+1+j]))*P_temps[j+1];
		}
		dec++;
		k--;
	}

	glm::vec3 test = P_temps[0];
	return P_temps[0];


	/*
	std::vector<glm::vec3> d;
	unsigned long p = _order;
	float x = u;
	std::vector<float> t = _knots;
	std::vector<glm::vec3> c = _controls;
	float alpha = 0.0f;

	for (unsigned long j = 0; j < p+1; ++j) {
		d.push_back(c[j+k-p]);	
	}
	

	for (unsigned long r = 1; r < p+1; ++r) {
		for (float j = p; j < -1; j += (r-1)) {
			alpha = (x-t[j+k-p])/(t[j+1+k-r]-t[j+k-p]);
			d[j] = (1.0f-alpha)*d[j-1]+alpha*d[j];
		}
	}

	glm::vec3 test = d[p];
	std::cout << test.x << " " << test.y << " " << test.z << std::endl;
	return test;
	*/


}

BSpline::~BSpline() {
	
}
