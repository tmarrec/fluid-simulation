#include "BSpline_tensor.h"

#include <iostream>

BSpline_tensor::BSpline_tensor(unsigned short order, std::vector<glm::vec3> controls, bool show_controls, float delta, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_window)
		: Shape(
			position,
			rotation,
			scale,
			geometry(order, controls, show_controls, delta),
			{RND(), RND(), RND()},
			{
				"shaders/vert.vert",
				"shaders/frag.frag"
			},
			main_window
		)
		, _order(order)
		, _controls(controls)
		, _delta(delta)
		, _show_controls(show_controls)
{
//	_knots = uniform_vector(controls.size()+order+1);
	set_geometry();
}

std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> BSpline_tensor::geometry(unsigned short order, std::vector<float> knots, std::vector<glm::vec3> controls, bool show_controls, float delta) {

	std::vector<GLfloat> vertex;
	std::vector<GLfloat> normals;
	std::vector<GLuint> indices;
	GLuint ind = 0;
	
	if (show_controls) {
		std::vector<float> flat;
		for (auto p : controls) {
			flat.push_back(p.x);
			flat.push_back(p.y);
			flat.push_back(p.z);
		}	
		vertex = flat;
		normals.assign(flat.size(), 1);
		indices.resize(controls.size());
		std::iota(std::begin(indices), std::end(indices), 0);
		ind = indices.size()-1;
	}

	for (float i = knots[order-1]; i < knots[controls.size()]; i += delta) {
		auto point = eval(i, order, knots, controls);
		vertex.push_back(point.x);
		vertex.push_back(point.y);
		vertex.push_back(point.z);
		normals.push_back(1);
		normals.push_back(1);
		normals.push_back(1);
		indices.push_back(ind++);
	}

	return {vertex, normals, indices};
}

std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> BSpline_tensor::geometry(unsigned short order, std::vector<glm::vec3> controls, bool show_controls, float delta) {
	std::vector<float> knots = uniform_vector(controls.size()+order+1);

	return geometry(order, knots, controls, show_controls, delta);
}

std::vector<float> BSpline_tensor::uniform_vector(unsigned short size) {
	std::vector<float> vec;
	for (unsigned short i = 0; i < size; ++i) {
		vec.push_back(i);
	}
	return vec;
}

glm::vec3 BSpline_tensor::eval(float u) {
//	return eval(u, _order, _knots, _controls);
}

glm::vec3 BSpline_tensor::eval(float u, unsigned short order, std::vector<float> knots, std::vector<glm::vec3> controls) {
	unsigned short dec = 0;
	unsigned short i = order;
	unsigned short m = order-1;

	while (u > knots[i]) {
		i++;
		dec++;
	}

	std::vector<glm::vec3> P_temps;
	unsigned short k = order;
	
	for(unsigned a = dec; a < dec + order; ++a)
		P_temps.emplace_back(controls[a]);

	for (unsigned short l = 0; l < m; ++l) {
		for (unsigned short j = 0; j < k-1; ++j) {
			P_temps[j] = ((knots[dec+k+j]-u)/(knots[dec+k+j]-knots[dec+1+j]))*P_temps[j]
						 +
						 ((u-knots[dec+1+j])/(knots[dec+k+j]-knots[dec+1+j]))*P_temps[j+1];
		}
		dec++;
		k--;
	}
	return P_temps[0];
}

BSpline_tensor::~BSpline_tensor() {
	
}

