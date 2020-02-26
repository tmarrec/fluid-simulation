#pragma once

#include "../Shape.h"

class Cube : public Shape {

public:
	Cube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
	~Cube();

	void draw(glm::vec3 view_position, glm::mat4 projection, float delta_time) override;


private:
	Shader _shader;

	glm::vec4 _color;
	void apply_color();
};
