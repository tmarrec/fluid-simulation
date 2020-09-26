#pragma once

#include "../Shape.h"
#include "../MainWindow.h"

class BSpline : public Shape {

public:
	BSpline(
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

	BSpline(
		unsigned short order,
		std::vector<glm::vec3> controls,
		bool show_controls,
		float delta,
		glm::vec3 position,
		glm::vec3 rotation,
		glm::vec3 scale,
		MainWindow * main_windown
	);

	~BSpline();

private:
	unsigned short _order; 
	std::vector<float> _knots;
	std::vector<glm::vec3> _controls; 
	float _delta;
	bool _show_controls;

	glm::vec3 eval(float u);
	glm::vec3 eval(float u, unsigned short order, std::vector<float> knots, std::vector<glm::vec3> controls);

	std::vector<float> uniform_vector(unsigned short size);

	std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> geometry(unsigned short order, std::vector<float> knots, std::vector<glm::vec3> controls, bool show_controls, float delta);
	std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> geometry(unsigned short order, std::vector<glm::vec3> controls, bool show_controls, float delta);

};
