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
			"Quad",
			{1.0f, 0.5f, 1.0f},
			{
				"shaders/vert.vert",
				"shaders/frag.frag"
			}
		)
{

}

Quad::~Quad() {
	
}

void Quad::draw(glm::vec3 view_position, glm::mat4 projection, float delta_time) {
	shader().use();
	apply_color();

	shader().set_mat4("model", get_model());
	shader().set_mat4("view", get_view(view_position));
	shader().set_mat4("projection", projection);

	shader().set_3f("_view_pos", view_position);
	shader().set_3f("_light_color", glm::vec3{1.0f, 1.0f, 1.0f});
	shader().set_3f("_light_pos", glm::vec3{0.0f, 0.0f, 10.0f});
	draw_vertex_quads();
}

void Quad::apply_color() {
	shader().set_3f("_object_color", color());
}

