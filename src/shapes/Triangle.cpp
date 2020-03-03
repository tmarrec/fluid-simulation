#include "Triangle.h"

#include <iostream>

Triangle::Triangle(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
		: Shape(
			position,
			rotation,
			scale,
			std::vector<GLfloat>{
				0.0f,		0.5f, 		0.0f, 		// Top Center
				0.5f,		-0.5f, 		0.0f, 		// Bottom Left
				-0.5f,		-0.5f, 		0.0f, 		// Bottom Right
			},
			std::vector<GLfloat>{
				0,		1,		0,
				1,		0,		1,
			},
			std::vector<GLuint>{
				0, 2, 1, // Triangle
			},
			"triangle",
			{1.0f, 0.5f, 1.0f}
		)
		, _shader {
			"../src/shaders/vert.vert",
			"../src/shaders/frag.frag"
		}
{

}

Triangle::~Triangle() {
	
}

void Triangle::draw(glm::vec3 view_position, glm::mat4 projection, float delta_time) {
	_shader.use();
	apply_color();

	rotate_test(delta_time);

	_shader.set_mat4("model", get_model());
	_shader.set_mat4("view", get_view(view_position));
	_shader.set_mat4("projection", projection);

	draw_vertex();
}

void Triangle::apply_color() {
	_shader.set_3f("_color", color());
}

