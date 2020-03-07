#pragma once


#include "../Shape.h"

class Model : public Shape {

public:
	Model(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
	~Model();

private:
	std::vector<GLfloat> get_vertices();
	std::vector<GLuint> get_indices();
	std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> get_geometry();
};
