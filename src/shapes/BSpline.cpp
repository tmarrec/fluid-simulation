#include "BSpline.h"

#include <iostream>

BSpline::BSpline(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_window)
		: Shape(
			position,
			rotation,
			scale,
			{
				{
					0.0f,		0.0f, 		0.0f, 		// Top Center
					1.0f,		1.0f, 		0.5f, 		// Bottom Left

					1.0f,		1.0f, 		0.5f, 		// Bottom Left
					2.0f,		1.0f, 		1.0f, 		// Bottom Right

				},
				{
					1,		1,		1,
					1,		1,		1,

					1,		1,		1,
					1,		1,		1,

				},
				{
					0, 1, 2,
					3
				},
			},
			{RND(), RND(), RND()},
			{
				"shaders/vert.vert",
				"shaders/frag.frag"
			},
			main_window
		)
{
	_order = 2;
	_controls = {
		{0, 0, 0},
		{1, 1, 0.5},
		{2, 1, 1}
	};
	_U = {0, 1, 2, 3, 4, 5}; // taille = _order + taille(controls) + 1

	set_geometry();

	for (unsigned short i = _U[_order-1]; i < _U[_controls.size()+1]; ++i) {
		eval(i);
	}
}

glm::vec3 BSpline::eval(float u) {
	unsigned short dec = 0;
	unsigned short i = _order;
	unsigned short k = _order;
	unsigned short m = _order-1;

	while (u > _U[i]) {
		i++;
		dec++;
	}
	std::vector<glm::vec3> P_temps;
	
	for(unsigned a = dec; a < dec + _order; ++a)
		P_temps.emplace_back(_controls[a]);

	for (unsigned short l = 0; l < m; ++l) {
		for (unsigned short j = 0; j < k-1; ++j) {
			P_temps[j] = (_U[dec+k+j]-u)/(_U[dec+k+j]-_U[dec+1+j])*P_temps[j]
						 +
						 (u-_U[dec+1+j])/(_U[dec+k+j]-_U[dec+1+j])*P_temps[j+1];
		}
		dec++;
		k--;
	}

	return P_temps[0];
}

BSpline::~BSpline() {
	
}

