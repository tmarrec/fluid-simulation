#pragma once

#include "../Shape.h"
#include "../MainWindow.h"
#include "BSpline.h"

class BSpline_tensor : public Shape {

public:
	BSpline_tensor(
		unsigned short order,
		std::vector<std::vector<glm::vec3>> controls,
		bool show_controls,
		float delta,
		glm::vec3 position,
		glm::vec3 rotation,
		glm::vec3 scale,
		MainWindow * main_windown
	);

	~BSpline_tensor();

private:
	std::vector<BSpline *> _leading_bsplines;
	std::vector<BSpline *> _generator_bsplines;
	float _delta;
	bool _show_controls;
	std::vector<std::vector<glm::vec3>> _controls;

	unsigned long _nb_gene;
	unsigned long _nb_point_gene;

	std::tuple<float, float> _range;

	void generate_leading_bsplines(unsigned short order);
	void generate_generator_bsplines(unsigned short order);
	std::vector<float> uniform_vector(unsigned short size);

	std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> geometry();

};
