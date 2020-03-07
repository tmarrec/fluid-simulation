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

void Cube::draw(glm::vec3 view_position, glm::mat4 projection, float delta_time,
				std::vector<std::shared_ptr<Entity>> lights) {
	shader().use();
	apply_color();

	//rotate_test(delta_time); 

	shader().set_mat4("model", get_model());
	shader().set_mat4("view", get_view(view_position));
	shader().set_mat4("projection", projection);

	shader().set_3f("_view_pos", view_position);
	shader().set_1i("_light_nb", lights.size());

	for (size_t i = 0; i < lights.size(); ++i) {
		auto light = std::static_pointer_cast<Light>(lights[i]);
		auto temp = std::string("_point_lights[") + std::to_string(i) + "].position";
		shader().set_3f(temp.c_str(), light->position());
		temp = std::string("_point_lights[") + std::to_string(i) + "].color";
		shader().set_3f(temp.c_str(), light->color());
	}

	draw_vertex();
}

void Cube::apply_color() {
	shader().set_3f("_object_color", color());
}

