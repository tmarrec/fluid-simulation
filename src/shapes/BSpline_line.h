#pragma once

#include "../Shape.h"
#include "../MainWindow.h"
#include "BSpline.h"

class BSpline_line : public Shape {

public:
	BSpline_line(
		unsigned short order,
		std::vector<float> knots,
		std::vector<glm::vec3> controls,
		bool show_controls,
		float delta,
		glm::vec3 position,
		glm::vec3 rotation,
		glm::vec3 scale,
		MainWindow * main_windown
	);

	BSpline_line(
		unsigned short order,
		std::vector<glm::vec3> controls,
		bool show_controls,
		float delta,
		glm::vec3 position,
		glm::vec3 rotation,
		glm::vec3 scale,
		MainWindow * main_windown
	);

	~BSpline_line();

private:
	BSpline * _bspline;
	float _delta;
	bool _show_controls;
	std::vector<glm::vec3> _controls;

	std::tuple<float, float> _range;

	std::vector<float> uniform_vector(unsigned short size);

	std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> geometry();

};
