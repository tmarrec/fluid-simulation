#pragma once

#include "../Shape.h"

class Quad : public Shape {

public:
	Quad(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
	~Quad();

private:
	std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> get_geometry();
};
