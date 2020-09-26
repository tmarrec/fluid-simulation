#include "BSpline.h"

#include <iostream>

BSpline::BSpline(unsigned short order, std::vector<float> knots, std::vector<glm::vec3> controls, bool show_controls, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_window)
		: Shape(
			position,
			rotation,
			scale,
			geometry(order, knots, controls, show_controls),
			{RND(), RND(), RND()},
			{
				"shaders/vert.vert",
				"shaders/frag.frag"
			},
			main_window
		)
		, _order(order)
		, _knots(knots)
		, _controls(controls)
		, _show_controls(show_controls)
{
	set_geometry();
}

std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> BSpline::geometry(unsigned short order, std::vector<float> knots, std::vector<glm::vec3> controls, bool show_controls) {

	/*
	for (unsigned short i = 0; i < _controls.size()+_order+1; ++i) {
		_knots.push_back(i);
	}
	*/

	std::vector<GLfloat> vertex;
	std::vector<GLfloat> normals;
	std::vector<GLuint> indices;
	GLuint ind = 0;
	
	if (show_controls) {
		std::vector<float> test;
		for (auto p : controls) {
			test.push_back(p.x);
			test.push_back(p.y);
			test.push_back(p.z);
		}	
		vertex = test;
		normals.assign(test.size(), 1);
		indices.resize(controls.size());
		std::iota(std::begin(indices), std::end(indices), 0);
	}

	/*
	for (float i = _knots[_order-1]; i < _knots[_controls.size()+1]; i += 0.1f) {
		auto point = eval(i);
		vertex.push_back(point.x);
		vertex.push_back(point.y);
		vertex.push_back(point.z);
		std::cout << glm::to_string(point) << std::endl;
		normals.push_back(1);
		normals.push_back(1);
		normals.push_back(1);
		indices.push_back(ind++);
	}
	*/
	/*
	return 
	{				{
					0.0f,		0.5f, 		0.0f, 		// Top Center
					0.5f,		-0.5f, 		0.0f, 		// Bottom Left
					-0.5f,		-0.5f, 		0.0f, 		// Bottom Right
				},
				{
					1,		1,		1,
					1,		1,		1,
					1,		1,		1,
				},
				{
					0, 2, 1, // Triangle
				},
	};
	*/
	return {vertex, normals, indices};
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
	
	for(unsigned a = dec; a < dec + _order; ++a)
		P_temps.emplace_back(_controls[a]);

	for (unsigned short l = 0; l < m; ++l) {
		for (unsigned short j = 0; j < k-1; ++j) {
			P_temps[j] = ((_knots[dec+k+j]-u)/(_knots[dec+k+j]-_knots[dec+1+j]))*P_temps[j]
						 +
						 ((u-_knots[dec+1+j])/(_knots[dec+k+j]-_knots[dec+1+j]))*P_temps[j+1];
		}
		dec++;
		k--;
	}
	std::cout << " ??? " << std::endl;
	return P_temps[0];
}

BSpline::~BSpline() {
	
}

