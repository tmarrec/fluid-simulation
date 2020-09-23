#include "Shape.h"
#include <iostream>

#include <chrono>

Shape::Shape(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
	std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> geometry, 
	glm::vec3 color, Shader shader, MainWindow * main_window)
	: Entity(position, rotation, scale, main_window)
	, _vertices{std::get<0>(geometry)}
	, _normals{std::get<1>(geometry)}
	, _indices{std::get<2>(geometry)}
	, _color{color}
	, _shader{shader}
{

}

Shape::~Shape(void) {
	glDeleteVertexArrays(1, &_VAO);
	glDeleteBuffers(1, &_VBO);
	glDeleteBuffers(1, &_NBO);
	glDeleteBuffers(1, &_EBO);
}

Entity_Type Shape::type() {
	return SHAPE;
}

void Shape::set_shader(Shader shader) {
	_shader = shader;
}

void Shape::set_geometry() {
	// Initialize the geometry
	glGenBuffers(1, &_VBO);
	glGenBuffers(1, &_NBO);
	glGenBuffers(1, &_EBO);

	glGenVertexArrays(1, &_VAO);

	glBindVertexArray(_VAO);
		
		glBindBuffer(GL_ARRAY_BUFFER, _VBO);
		glBufferData(
			GL_ARRAY_BUFFER,
			_vertices.size()*sizeof(GLfloat),
			_vertices.data(),
			GL_STATIC_DRAW
		);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, _NBO);
		glBufferData(
			GL_ARRAY_BUFFER,
			_normals.size()*sizeof(GLfloat),
			_normals.data(),
			GL_STATIC_DRAW
		);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _EBO);
		glBufferData(
			GL_ELEMENT_ARRAY_BUFFER,
			_indices.size()*sizeof(GLuint),
			_indices.data(),
			GL_STATIC_DRAW
		);

	glBindVertexArray(0);
}

std::vector<GLuint> Shape::indices() const {
	return _indices;
}

const glm::vec3 Shape::color() const {
	return _color;
}

void Shape::draw_vertex() {
	glBindVertexArray(_VAO);
	glDrawElements(GL_LINES, indices().size(), GL_UNSIGNED_INT, nullptr);
	glBindVertexArray(0);
}

// Utilise le shader et envoie les uniforms 
void Shape::use_shader(glm::mat4 view, glm::mat4 projection,
				std::vector<std::shared_ptr<Entity>> lights) {
	shader().use();
	shader().set_3f("_object_color", color());

	shader().set_mat4("model", get_model());
	shader().set_mat4("view", view);
	shader().set_mat4("projection", projection);
	
	uint64_t time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	float time_sin = sin(time/50000);
	shader().set_1f("time", time_sin);

	shader().set_3f("_view_pos", view[3]);
	shader().set_1i("_light_nb", lights.size());

	// Envoie les uniforms pour toutes les lumieres
	for (size_t i = 0; i < lights.size(); ++i) {
		auto light = std::static_pointer_cast<Light>(lights[i]);
		auto temp = std::string("_point_lights[") + std::to_string(i) + "].position";
		shader().set_3f(temp.c_str(), light->position());
		temp = std::string("_point_lights[") + std::to_string(i) + "].color";
		shader().set_3f(temp.c_str(), light->color());
		temp = std::string("_point_lights[") + std::to_string(i) + "].intensity";
		shader().set_1f(temp.c_str(), light->intensity());
		temp = std::string("_point_lights[") + std::to_string(i) + "].constant";
		shader().set_1f(temp.c_str(), 1.0f);
		temp = std::string("_point_lights[") + std::to_string(i) + "].linear";
		shader().set_1f(temp.c_str(), 0.0014f);
		temp = std::string("_point_lights[") + std::to_string(i) + "].quadratic";
		shader().set_1f(temp.c_str(), 0.000007f);
	}
}

void Shape::draw(glm::mat4 view, glm::mat4 projection, float delta_time,
				std::vector<std::shared_ptr<Entity>> lights) {
	//rotate_test(delta_time); 
	std::ignore = delta_time;
	use_shader(view, projection, lights);
	draw_vertex();
}

Shader & Shape::shader() {
	return _shader;
}

