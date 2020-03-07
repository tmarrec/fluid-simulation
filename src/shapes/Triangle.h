#pragma once

#include "../Shape.h"

class Triangle : public Shape {

public:
	Triangle(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
	~Triangle();

private:
	std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> get_geometry();
};
