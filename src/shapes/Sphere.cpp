#include "Sphere.h"

#include <iostream>

Sphere::Sphere(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, glm::vec2 faces)
		: Shape(
			position,
			rotation,
			scale,
			generate_vertices(faces),
			generate_vertices(faces),
			generate_indices(faces),
			"Sphere",
			{0.0f, 0.0f, 1.0f},
			{
				"shaders/error.vert",
				"shaders/frag.frag"
			}	
		)
{

}

std::vector<GLfloat> Sphere::generate_vertices(glm::vec2 faces) {
	std::vector<GLfloat> vertices;
	float x, y, z, xy;

	unsigned int longitude_count = faces[0];
	unsigned int latitude_count = faces[1];

	float longitude_step = 2*glm::pi<float>() / longitude_count;
	float latitude_step = glm::pi<float>() / latitude_count;
	float longitude_angle, latitude_angle;

	for (unsigned int i = 0; i <= latitude_count; ++i) {
		latitude_angle = glm::pi<float>() / 2 - i * latitude_step;
		xy = cosf(latitude_angle); 
		z = sinf(latitude_angle); 
		for (unsigned int j = 0; j <= longitude_count; ++j) {
			longitude_angle = j*longitude_step;

			x = xy * cosf(longitude_angle);
			y = xy * sinf(longitude_angle);

			vertices.push_back(x);
			vertices.push_back(y);
			vertices.push_back(z);
		}
	}
	return vertices;
}

std::vector<GLuint> Sphere::generate_indices(glm::vec2 faces) {
	std::vector<GLuint> indices;
	int k1, k2;
	unsigned int longitude_count = faces[0];
	unsigned int latitude_count = faces[1];

	for (unsigned int i = 0; i < latitude_count; ++i) {
		k1 = i * (longitude_count+1);
		k2 = k1 + longitude_count + 1;

		for (unsigned int j = 0; j < longitude_count; ++j, ++k1, ++k2) {
			if (i != 0) {
				indices.push_back(k1);
				indices.push_back(k2);
				indices.push_back(k1+1);
			}
			if (i != latitude_count-1) {
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

void Sphere::apply_color() {
	shader().set_3f("_object_color", color());
}

