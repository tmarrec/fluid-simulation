#pragma once

#include <vector>
#include <glm/glm.hpp>

class BSpline {

public:
	BSpline(
		unsigned short order,
		std::vector<float> knots,
		std::vector<glm::vec3> controls
	);

	glm::vec3 eval(float u);

	~BSpline();

private:
	unsigned short _order;
	std::vector<float> _knots;
	std::vector<glm::vec3> _controls;
};
