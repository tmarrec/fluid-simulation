#pragma once

#include "../Shape.h"

class Triangle : public Shape {

public:
	Triangle(glm::vec3 position, glm::vec3 rotation, glm::vec2 scale);
	~Triangle();

	virtual void draw(glm::vec3 view_position, glm::mat4 projection, float delta_time);


private:
	Shader _shader;

	glm::vec4 _color;
	void apply_color();
};
