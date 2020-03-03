#include "Cube.h"


#include <iostream>

Cube::Cube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
		: Shape(
			position,
			rotation,
			scale,
			std::vector<GLfloat>{
				// front
				-0.5f,		0.5f, 		-0.5f, 		// Top Left
				0.5f,		0.5f, 		-0.5f, 		// Top Right
				-0.5f,		-0.5f, 		-0.5f, 		// Bottom Left
				0.5f,		-0.5f, 		-0.5f, 		// Bottom Right

				// back
				-0.5f,		0.5f, 		0.5f, 		// Top Left
				0.5f,		0.5f, 		0.5f, 		// Top Right
				-0.5f,		-0.5f, 		0.5f, 		// Bottom Left
				0.5f,		-0.5f, 		0.5f, 		// Bottom Right
			},
			std::vector<GLfloat>{
				1,		0,		1,
				1,		0,		1,
				1,		0,		1,
				0,		1,		0,
				1,		0,		1,
				1,		0,		1,
				1,		0,		1,
				1,		0,		1,
			},
			std::vector<GLuint>{
				0, 2, 3,
				0, 3, 1,
				1, 3, 7,
				1, 7, 5,
				5, 7, 6,
				5, 6, 4,
				4, 6, 2,
				4, 2, 0,
				6, 2, 3,
				6, 3, 7,
				4, 0, 1,
				4, 1, 5,
			},
			"cube",
			{1.0f, 1.0f, 1.0f}
		)
		, _shader {
			"../src/shaders/vert.vert",
			"../src/shaders/frag.frag"
		}
{

}

Cube::~Cube() {
	
}

void Cube::draw(glm::vec3 view_position, glm::mat4 projection, float delta_time) {
	_shader.use();
	apply_color();

	rotate_test(delta_time);

	_shader.set_mat4("model", get_model());
	_shader.set_mat4("view", get_view(view_position));
	_shader.set_mat4("projection", projection);

	draw_vertex();
}

void Cube::apply_color() {
	_shader.set_3f("_color", color());
}

