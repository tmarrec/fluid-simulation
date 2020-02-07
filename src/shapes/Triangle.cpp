#include "Triangle.h"

#include <iostream>
#include <string>
#include <fstream>
#include <streambuf>

#include <QDateTime>

Triangle::Triangle(glm::vec3 position, glm::vec3 rotation, glm::vec2 scale)
		: Entity(position, rotation, scale)
		, _shader {
			"../src/shaders/triangle.vert",
			"../src/shaders/triangle.frag"
		}
		, _color {1.0f, 0.5f, 1.0f, 1.0f}
		, _position {position}
		, _rotation {rotation}
		, _scale {scale}
{
	_vertices = {
		0.0f,		0.5f, 		0.0f, 		// Top Center
		0.5f,		-0.5f, 		0.0f, 		// Bottom Left
		-0.5f,		-0.5f, 		0.0f, 		// Bottom Right
	};

	_normals = {
		0,		1,		0,
		1,		0,		1,
	};

	_indices = {
		0, 2, 1, // Triangle
	};

	// Initialize the geometry
	// 1. Generate geometry buffers
	glGenBuffers(1, &_vbo);
	glGenBuffers(1, &_nbo);
	glGenBuffers(1, &_ebo);
	glGenVertexArrays(1, &_vao);

	// 2. Bind Vertex Array Object
	glBindVertexArray(_vao);
		
		// 3. Copy our vertices array in a buffer for OpenGL to use
		glBindBuffer(GL_ARRAY_BUFFER, _vbo);
		glBufferData(
			GL_ARRAY_BUFFER,
			_vertices.size()*sizeof(GLfloat),
			_vertices.data(),
			GL_STATIC_DRAW
		);

		// 4. Then set our vertex attributes pointers
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(0);

		// 5. Copy our normals array in a buffer for OpenGL to use
		glBindBuffer(GL_ARRAY_BUFFER, _nbo);
		glBufferData(
			GL_ARRAY_BUFFER,
			_normals.size()*sizeof(GLfloat),
			_normals.data(),
			GL_STATIC_DRAW
		);

		// 6. Copy our vertices array in a buffer for OpenGL to use
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(1);

		// 7. Copy our index array in a element buffer for OpenGL to use
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _ebo);
		glBufferData(
			GL_ELEMENT_ARRAY_BUFFER,
			_indices.size()*sizeof(GLfloat),
			_indices.data(),
			GL_STATIC_DRAW
		);

	// 8. Unbind the VAO
	glBindVertexArray(0);
}

Triangle::~Triangle() {
	
}

void Triangle::draw(glm::vec3 view_position, glm::mat4 projection, float delta_time) {
	_shader.use();

	apply_color();

	_rotation.x += 150 * delta_time;

	// Matrice model pour definir la position
	// et la rotation de l'objet dans l'espace
	glm::mat4 model {1.0f};
	model = glm::translate(model, _position);
	model = glm::rotate(model, glm::radians(_rotation.x), glm::vec3(1.0f, 0.0f, 0.0f));
	model = glm::rotate(model, glm::radians(_rotation.y), glm::vec3(0.0f, 1.0f, 0.0f));
	model = glm::rotate(model, glm::radians(_rotation.z), glm::vec3(0.0f, 0.0f, 1.0f));
	model = glm::scale(model, glm::vec3{_scale, 1.0f});

	glm::mat4 view {1.0f};
	view = glm::translate(view, view_position);

	_shader.set_mat4("model", model);
	_shader.set_mat4("view", view);
	_shader.set_mat4("projection", projection);

	glBindVertexArray(_vao);
	glDrawElements(GL_TRIANGLES, _indices.size(), GL_UNSIGNED_INT, nullptr);
	glBindVertexArray(0);
}

void Triangle::apply_color() {
	_shader.set_4f("_color", _color);
}

void Triangle::set_position(glm::vec3 position) {
	_position = position;
}

void Triangle::set_rotation(glm::vec3 rotation) {
	_rotation = rotation;
}

glm::vec3 Triangle::position() const {
	return _position;
}

glm::vec3 Triangle::rotation() const {
	return _rotation;
}
