#include "Sphere.h"

#include <iostream>

Sphere::Sphere(glm::vec3 position, glm::vec3 rotation, glm::vec2 scale)
		: Shape(
			position,
			rotation,
			scale,
			generate_vertices(),
			std::vector<GLfloat>{
				0,		1,		0,
				1,		0,		1,
			},
			generate_indices()	
		)
		, _shader {
			"../src/shaders/triangle.vert",
			"../src/shaders/triangle.frag"
		}
		, _color {0.0f, 0.0f, 1.0f, 1.0f}
{

}

std::vector<GLfloat> Sphere::generate_vertices() {
	std::vector<GLfloat> vertices;
	float x, y, z, xy;
	float nx, ny, nz;

	// temp
	float sector_count = 10;
	float stack_count = 10;
	float radius = 1;


	float sector_step = 2*glm::pi<float>() / sector_count;
	float stack_step = glm::pi<float>() / stack_count;
	float sector_angle, stack_angle;

	for (int i = 0; i <= stack_count; ++i) {
		stack_angle = glm::pi<float>() / 2 - i * stack_step;
		xy = radius * cosf(stack_angle); 
		z = radius * sinf(stack_angle); 
		for (int j = 0; j <= sector_count; ++j) {
			sector_angle = j*sector_step;

			x = xy * cosf(sector_angle);
			y = xy * sinf(sector_angle);

			std::cout << x << " " << y << " " << z << std::endl;
			vertices.push_back(x);
			vertices.push_back(y);
			vertices.push_back(z);
		}
	}
	return vertices;
}

std::vector<GLuint> Sphere::generate_indices() {
	std::vector<GLuint> indices;
	int k1, k2;

	// temp
	float sector_count = 10;
	float stack_count = 10;
	float radius = 1;

	for (int i = 0; i < stack_count; ++i) {
		k1 = i * (sector_count+1);
		k2 = k1 + sector_count + 1;

		for (int j = 0; j < sector_count; ++j, ++k1, ++k2) {
			if (i != 0) {
				indices.push_back(k1);
				indices.push_back(k2);
				indices.push_back(k1+1);
			}
			if (i != stack_count-1) {
				indices.push_back(k1+1);
				indices.push_back(k2);
				indices.push_back(k2+1);
			}
		}
	}
	return indices;
}

Sphere::~Sphere() {
	
}

void Sphere::draw(glm::vec3 view_position, glm::mat4 projection, float delta_time) {
	_shader.use();
	apply_color();

	rotate_test(delta_time);

	_shader.set_mat4("model", get_model());
	_shader.set_mat4("view", get_view(view_position));
	_shader.set_mat4("projection", projection);

	draw_vertex();
}

void Sphere::apply_color() {
	_shader.set_4f("_color", _color);
}

