#include "Triangle.h"
#include "../OpenGL.h"
#include "../Shader.h"

#include <iostream>
#include <string>
#include <fstream>
#include <streambuf>

#include <QDateTime>

Triangle::Triangle(int w, int h)
		: OpenGL {w, h}
		, _shader {
			"../src/shaders/triangle.vert",
			"../src/shaders/triangle.frag"
		}
		, _color {1.0f, 0.5f, 1.0f, 1.0f}
{
	std::cout << "Triangle construct" << std::endl;
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

	_shader.use();
}

Triangle::~Triangle() {
	
}

void Triangle::draw() {
	OpenGL::draw();

	_shader.use();

	apply_color();

	float time = sin(QDateTime::currentMSecsSinceEpoch()) / 2.0f + 0.5f;
	
	test += 100 * delta_time();

	glm::mat4 _model {1.0f};
	_model = glm::rotate(_model, glm::radians(test), glm::vec3(0.0f, 1.0f, 0.0f));

	glm::mat4 _view {1.0f};
	_view = glm::translate(_view, glm::vec3(0.0f, 0.0f, -2.0f));

	_shader.set_mat4("model", _model);
	_shader.set_mat4("view", _view);
	_shader.set_mat4("projection", projection());

	glBindVertexArray(_vao);
	glDrawElements(GL_TRIANGLES, _indices.size(), GL_UNSIGNED_INT, nullptr);
	glBindVertexArray(0);
}

void Triangle::apply_color() {
	_shader.set_4f("_color", _color);
}
