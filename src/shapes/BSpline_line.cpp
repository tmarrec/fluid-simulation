#include "BSpline_line.h"

#include <iostream>

BSpline_line::BSpline_line(unsigned short order, std::vector<float> knots, std::vector<glm::vec3> controls, bool show_controls, float delta, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_window)
		: Shape(
			position,
			rotation,
			scale,
			{},
			{RND(), RND(), RND()},
			{
				"shaders/vert.vert",
				"shaders/frag.frag"
			},
			main_window
		)
		, _delta(delta)
		, _show_controls(show_controls)
		, _controls(controls)
		, _range(knots[order-1], knots[controls.size()])
{
	_bspline = new BSpline(order, knots, controls);
	set_vert_norm_indi(geometry());
	set_geometry();
}

BSpline_line::BSpline_line(unsigned short order, std::vector<glm::vec3> controls, bool show_controls, float delta, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_window)
		: Shape(
			position,
			rotation,
			scale,
			{},
			{RND(), RND(), RND()},
			{
				"shaders/vert.vert",
				"shaders/frag.frag"
			},
			main_window
		)
		, _delta(delta)
		, _show_controls(show_controls)
		, _controls(controls)
{
	auto knots = uniform_vector(controls.size()+order+1);
	_range = std::tuple<float,float>(knots[order-1], knots[controls.size()]);
	_bspline = new BSpline(order, knots, controls);
	set_vert_norm_indi(geometry());
	set_geometry();
}

std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> BSpline_line::geometry() {

	std::vector<GLfloat> vertex;
	std::vector<GLfloat> normals;
	std::vector<GLuint> indices;
	GLuint ind = 0;
	
	if (_show_controls) {
		std::vector<float> flat;
		for (auto p : _controls) {
			flat.push_back(p.x);
			flat.push_back(p.y);
			flat.push_back(p.z);
		}	
		vertex = flat;
		normals.assign(flat.size(), 1);
		indices.resize(_controls.size());
		std::iota(std::begin(indices), std::end(indices), 0);
		ind = indices.size()-1;
	}

	for (float i = std::get<0>(_range); i < std::get<1>(_range); i += _delta) {
		auto point = _bspline->eval(i);
		std::cout << glm::to_string(point) << std::endl;
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

std::vector<float> BSpline_line::uniform_vector(unsigned short size) {
	std::vector<float> vec;
	for (unsigned short i = 0; i < size; ++i) {
		vec.push_back(i);
	}
	return vec;
}

BSpline_line::~BSpline_line() {
	
}

