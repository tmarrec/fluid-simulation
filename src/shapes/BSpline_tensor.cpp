#include "BSpline_tensor.h"

#include <iostream>

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
	generate_leading_bsplines(order);
	generate_generator_bsplines(order);

	_nb_gene = _generator_bsplines.size();
	// TODO dernier probleme
	std::cout << "a: " << std::get<0>(_generator_bsplines.front()->range()) << std::endl;
	std::cout << "b: " << std::get<1>(_generator_bsplines.front()->range()) << std::endl;
	std::cout << "c: " << _delta << std::endl;
	_nb_point_gene = std::floor((std::get<1>(_generator_bsplines.front()->range())-std::get<0>(_generator_bsplines.front()->range()))/_delta);
	if (_nb_point_gene%2 == 1) {
		_nb_point_gene++;
	}
	std::cout << "res: " << _nb_point_gene << std::endl;

	std::cout << "LE DEBUG" << std::endl;

	set_vert_norm_indi(geometry());
	set_geometry();
}

void BSpline_tensor::generate_leading_bsplines(unsigned short order) {
	for (auto bs : _controls) {
		auto knots = uniform_vector(bs.size()+order+1);
		_leading_bsplines.push_back(new BSpline(order, knots, bs));
	}
}

void BSpline_tensor::generate_generator_bsplines(unsigned short order) {
	for (float i = std::get<0>(_leading_bsplines.front()->range()); i < std::get<1>(_leading_bsplines.front()->range()); i += _delta) {
		std::vector<glm::vec3> generator_controls;
		for (auto lbs : _leading_bsplines) {
			auto point = lbs->eval(i);
			generator_controls.push_back(point);
		}
		auto knots = uniform_vector(generator_controls.size()+order+1);
		auto gbs = new BSpline(order, knots, generator_controls);
		_generator_bsplines.push_back(gbs);
	}
}

std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> BSpline_tensor::geometry() {

	std::vector<GLfloat> vertex;
	std::vector<GLfloat> normals;
	std::vector<GLuint> indices;
	
	/* TODO fix ca pour montrer tout les controls des bsplines
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
	*/

	for (auto gbs : _generator_bsplines) {
		for (float i = std::get<0>(_generator_bsplines.front()->range()); i < std::get<1>(_generator_bsplines.front()->range()); i += _delta) {
			auto point = gbs->eval(i);

			vertex.push_back(point.x);
			vertex.push_back(point.y);
			vertex.push_back(point.z);
			normals.push_back(1);
			normals.push_back(1);
			normals.push_back(1);
		}
	}


	std::cout << "_nb_gene: " << _nb_gene << std::endl;
	std::cout << "_nb_point_gene: " << _nb_point_gene << std::endl;

	for (unsigned long i = 0; i < _nb_gene-1; ++i) {
		for (unsigned long j = 0; j < _nb_point_gene-1; ++j) {
			unsigned long p = j+i*_nb_point_gene;
			// First face triangle
			indices.push_back(p);
			indices.push_back(p+_nb_point_gene);
			indices.push_back(p+1);

			// Second face triangle
			indices.push_back(p+_nb_point_gene);
			indices.push_back(p+_nb_point_gene+1);
			indices.push_back(p+1);
		}
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

