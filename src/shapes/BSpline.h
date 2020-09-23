#pragma once

#include "../Shape.h"
#include "../MainWindow.h"

class BSpline : public Shape {

public:
	BSpline(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_windown);
	~BSpline();

private:
	unsigned short _order; 
	std::vector<float> _U;
	std::vector<glm::vec3> _controls; 

	glm::vec3 eval(float u);

};
