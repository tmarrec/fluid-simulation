#include "Cube.h"


#include <iostream>

Cube::Cube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
		: Shape(
			position,
			rotation,
			scale,
			std::vector<GLfloat>{
				-0.5f, -0.5f, 0.5f,
				0.5f, -0.5f, 0.5f,
				0.5f, 0.5f, 0.5f,
				-0.5f, 0.5f, 0.5f,

				0.5f, 0.5f, 0.5f,
				0.5f, 0.5f, -0.5f,
				0.5f, -0.5f, -0.5f,
				0.5f, -0.5f, 0.5f,

				-0.5f, -0.5f, -0.5f,
				0.5f, -0.5f, -0.5f,
				0.5f, 0.5f, -0.5f,
				-0.5f, 0.5f, -0.5f,

				-0.5f, -0.5f, -0.5f,
				-0.5f, -0.5f, 0.5f,
				-0.5f, 0.5f, 0.5f,
				-0.5f, 0.5f, -0.5f,
				
				0.5f, 0.5f, 0.5f,
				-0.5f, 0.5f, 0.5f,
				-0.5f, 0.5f, -0.5f,
				0.5f, 0.5f, -0.5f,

				-0.5f, -0.5f, -0.5f,
				0.5f, -0.5f, -0.5f,
				0.5f, -0.5f, 0.5f,
				-0.5f, -0.5f, 0.5f,
			},
			std::vector<GLfloat>{
				0.0f, 0.0f, 1.0f,
				0.0f, 0.0f, 1.0f,
				0.0f, 0.0f, 1.0f,
				0.0f, 0.0f, 1.0f,

				1.0f, 0.0f, 0.0f,
				1.0f, 0.0f, 0.0f,
				1.0f, 0.0f, 0.0f,
				1.0f, 0.0f, 0.0f,

				0.0f, 0.0f, -1.0f,
				0.0f, 0.0f, -1.0f,
				0.0f, 0.0f, -1.0f,
				0.0f, 0.0f, -1.0f,

				-1.0f, 0.0f, 0.0f,
				-1.0f, 0.0f, 0.0f,
				-1.0f, 0.0f, 0.0f,
				-1.0f, 0.0f, 0.0f,

				0.0f, 1.0f, 0.0f,
				0.0f, 1.0f, 0.0f,
				0.0f, 1.0f, 0.0f,
				0.0f, 1.0f, 0.0f,

				0.0f, -1.0f, 0.0f,
				0.0f, -1.0f, 0.0f,
				0.0f, -1.0f, 0.0f,
				0.0f, -1.0f, 0.0f,

			},
			std::vector<GLuint>{
				0,  1,  2,  0,  2,  3,   //front
 				4,  5,  6,  4,  6,  7,   //right
 				8,  9,  10, 8,  10, 11,  //back
 				12, 13, 14, 12, 14, 15,  //left
 				16, 17, 18, 16, 18, 19,  //upper
 				20, 21, 22, 20, 22, 23 	 //bottom
			},
			"Cube",
			{1.0f, 1.0f, 0.0f},
			{ // TODO c'est ignoble
				"shaders/vert.vert",
				"shaders/frag.frag"
			}
		)
{

}

Cube::~Cube() {
	
}

void Cube::draw(glm::vec3 view_position, glm::mat4 projection, float delta_time) {
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

void Cube::apply_color() {
	shader().set_3f("_object_color", color());
}

