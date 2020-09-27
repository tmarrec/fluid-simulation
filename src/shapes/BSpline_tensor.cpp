#include "BSpline_tensor.h"

#include <iostream>

BSpline_tensor::BSpline_tensor(unsigned short order, std::vector<float> knots, std::vector<std::vector<glm::vec3>> controls, bool show_controls, float delta, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_window)
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
		, _range(knots[order-1]+delta, knots[controls.size()-delta])
{
	generate_leading_bsplines(controls, order, knots);
	generate_generator_bsplines(order, knots);
	set_vert_norm_indi(geometry());
	set_geometry();
}

BSpline_tensor::BSpline_tensor(unsigned short order, std::vector<std::vector<glm::vec3>> controls, bool show_controls, float delta, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_window)
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
	_range = std::tuple<float,float>(knots[order-1]+delta, knots[controls.size()-_delta]);
	generate_leading_bsplines(controls, order, knots);
	generate_generator_bsplines(order, knots);
	set_vert_norm_indi(geometry());
	set_geometry();
}

void BSpline_tensor::generate_leading_bsplines(std::vector<std::vector<glm::vec3>> controls, unsigned short order, std::vector<float> knots) {
	for (auto bs : controls) {
		_leading_bsplines.push_back(new BSpline(order, knots, bs));
	}
}

void BSpline_tensor::generate_generator_bsplines(unsigned short order, std::vector<float> knots) {
	for (float i = std::get<0>(_range); i < std::get<1>(_range); i += _delta) {
		std::vector<glm::vec3> generator_controls;
		for (auto lbs : _leading_bsplines) {
			auto point = lbs->eval(i);
			generator_controls.push_back(point);
			std::cout << glm::to_string(point) << std::endl;
		}
		std::cout << std::endl;
		auto gbs = new BSpline(order, knots, generator_controls);
		_generator_bsplines.push_back(gbs);
	}
	std::cout << "################################" << std::endl;
}

std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> BSpline_tensor::geometry() {

	std::vector<GLfloat> vertex;
	std::vector<GLfloat> normals;
	std::vector<GLuint> indices;
	GLuint ind = 0;
	
	/* TODO fix ca pour montrer tout les controls des bsplines
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
	*/

	unsigned long nb_point_line = 0;
	for (auto gbs : _generator_bsplines) {
		nb_point_line = 0;
		for (float i = std::get<0>(_range); i < std::get<1>(_range); i += _delta) {
			nb_point_line++;
			auto point = gbs->eval(i);

			std::cout << glm::to_string(point) << std::endl;
			vertex.push_back(point.x);
			vertex.push_back(point.y);
			vertex.push_back(point.z);
			normals.push_back(1);
			normals.push_back(1);
			normals.push_back(1);
		}
		std::cout << std::endl;
	}
	unsigned long nb_gene = _generator_bsplines.size();

	std::cout << "nb_gene: " << nb_gene << std::endl;

	//unsigned long nb_point_line = std::floor((std::get<1>(_range)-std::get<0>(_range))/_delta);
			
	std::cout << "nb_point_line: " << nb_point_line << std::endl;

	for (unsigned long i = 0; i < nb_gene-1; ++i) {
		for (unsigned long j = 0; j < nb_point_line-1; ++j) {
			unsigned long p = j+i*nb_gene;
			// First face triangle
			indices.push_back(p);
			std::cout << indices.back() << " ";
			indices.push_back(p+nb_point_line);
			std::cout << indices.back() << " ";
			indices.push_back(p+1);
			std::cout << indices.back() << " ";

			std::cout << std::endl;
			// Second face triangle
			indices.push_back(p+nb_point_line);
			std::cout << indices.back() << " ";
			indices.push_back(p+nb_point_line+1);
			std::cout << indices.back() << " ";
			indices.push_back(p+1);
			std::cout << indices.back() << " ";
			std::cout << std::endl;
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	return {vertex, normals, indices};
}

std::vector<float> BSpline_tensor::uniform_vector(unsigned short size) {
	std::vector<float> vec;
	for (unsigned short i = 0; i < size; ++i) {
		vec.push_back(i);
	}
	return vec;
}

BSpline_tensor::~BSpline_tensor() {
	
}

