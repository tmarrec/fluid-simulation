#include "Quad.h"

#include <iostream>

Quad::Quad(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
		: Shape(
			position,
			rotation,
			scale,
			std::vector<GLfloat>{
				0.0f,		0.0f, 		0.0f, 		// Top Left
				0.5f,		-0.5f, 		0.0f, 		// Bottom Left
				-0.5f,		-0.5f, 		0.0f, 		// Bottom Right
				0.0f,		-0.5f, 		0.0f, 		// Top Right
			},
			std::vector<GLfloat>{
				0,		1,		0,
				1,		0,		1,
			},
			std::vector<GLuint>{
				0, 1, 2, 3 // Quad
			},
			"quad"
		)
		, _shader {
			"../src/shaders/triangle.vert",
			"../src/shaders/triangle.frag"
		}
		, _color {1.0f, 0.5f, 1.0f, 1.0f}
{

}

Quad::~Quad() {
	
}

void Quad::draw(glm::vec3 view_position, glm::mat4 projection, float delta_time) {
	_shader.use();
	apply_color();

	_shader.set_mat4("model", get_model());
	_shader.set_mat4("view", get_view(view_position));
	_shader.set_mat4("projection", projection);

	draw_vertex_quads();
}

void Quad::apply_color() {
	_shader.set_4f("_color", _color);
}

