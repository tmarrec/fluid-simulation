#pragma once

#include "../Shape.h"

class Sphere : public Shape {

public:
	Sphere(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, glm::vec2 faces);
	~Sphere();

private:
	std::vector<GLfloat> generate_vertices(glm::vec2 faces);
	std::vector<GLuint> generate_indices(glm::vec2 faces);
};
