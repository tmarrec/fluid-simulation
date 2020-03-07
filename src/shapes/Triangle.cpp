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
				1,		1,		1,
				1,		1,		1,
				1,		1,		1,
			},
			std::vector<GLuint>{
				0, 2, 1, // Triangle
			},
			"Triangle",
			{1.0f, 0.5f, 1.0f},
			{
				"shaders/vert.vert",
				"shaders/frag.frag"
			}
		)
{

}

Triangle::~Triangle() {
	
}

void Triangle::draw(glm::vec3 view_position, glm::mat4 projection, float delta_time,
				std::vector<std::shared_ptr<Entity>> lights) {
	shader().use();
	apply_color();

	rotate_test(delta_time);

	shader().set_mat4("model", get_model());
	shader().set_mat4("view", get_view(view_position));
	shader().set_mat4("projection", projection);

	shader().set_3f("_view_pos", view_position);
	shader().set_3f("_light_color", glm::vec3{1.0f, 1.0f, 1.0f});
	shader().set_3f("_light_pos", glm::vec3{0.0f, 0.0f, 10.0f});

	draw_vertex();
}

void Triangle::apply_color() {
	shader().set_3f("_object_color", color());
}

