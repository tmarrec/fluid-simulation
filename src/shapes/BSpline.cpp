#include "BSpline.h"

#include <iostream>

BSpline::BSpline(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_window)
		: Shape(
			position,
			rotation,
			scale,
			geometry(),
			{RND(), RND(), RND()},
			{
				"shaders/vert.vert",
				"shaders/frag.frag"
			},
			main_window
		)
{
	set_geometry();
}

std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> BSpline::geometry() {
		
	_order = 3;
	_controls = {
		{0, 0, 0},
		{1, 1, 0},
		{2, -1, 0},
		{3, 0, 0},
		{4, 1, 0},
		{5, -1, 0},
		{6, 0, 0},
	};

	for (unsigned short i = 0; i < _controls.size()+_order+1; ++i) {
		_U.push_back(i);
	}

	std::vector<GLfloat> vertex;
	std::vector<GLfloat> normals;
	std::vector<GLuint> indices;
	GLuint ind = 0;
	for (float i = _U[_order-1]; i < _U[_controls.size()+1]; i += 0.1f) {
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

	return {vertex, normals, indices};
}

glm::vec3 BSpline::eval(float u) {
	unsigned short dec = 0;
	unsigned short i = _order;
	unsigned short m = _order-1;

	while (u > _U[i]) {
		i++;
		dec++;
	}

	std::vector<glm::vec3> P_temps;
	unsigned short k = _order;
	
	for(unsigned a = dec; a < dec + _order; ++a)
		P_temps.emplace_back(_controls[a]);

	for (unsigned short l = 0; l < m; ++l) {
		for (unsigned short j = 0; j < k-1; ++j) {
			P_temps[j] = ((_U[dec+k+j]-u)/(_U[dec+k+j]-_U[dec+1+j]))*P_temps[j]
						 +
						 ((u-_U[dec+1+j])/(_U[dec+k+j]-_U[dec+1+j]))*P_temps[j+1];
		}
		dec++;
		k--;
	}
	std::cout << " ??? " << std::endl;
	return P_temps[0];
}

BSpline::~BSpline() {
	
}

