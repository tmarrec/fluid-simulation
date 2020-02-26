#pragma once

#include "../Shape.h"

class Sphere : public Shape {

public:
	Sphere(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, glm::vec2 faces);
	~Sphere();

	void draw(glm::vec3 view_position, glm::mat4 projection, float delta_time) override;


private:
	Shader _shader;

	glm::vec4 _color;
	void apply_color();
	std::vector<GLfloat> generate_vertices(glm::vec2 faces);
	std::vector<GLuint> generate_indices(glm::vec2 faces);
};
