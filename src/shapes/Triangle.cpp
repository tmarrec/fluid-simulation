#include "Triangle.h"
#include "../OpenGL.h"
#include "../Shader.h"

#include <iostream>
#include <string>
#include <fstream>
#include <streambuf>

#include <QDateTime>

Triangle::Triangle(int w, int h)
		: OpenGL(w, h)
		, shader {
			"../src/shaders/triangle.vert",
			"../src/shaders/triangle.frag"
		}
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

	shader.use();


}

Triangle::~Triangle() {
	
}

void Triangle::draw() {
	OpenGL::draw();

	
	//glUseProgram(_program);
	
	float time = sin(QDateTime::currentMSecsSinceEpoch()) / 2.0f + 0.5f;
	std::cout << time << std::endl;
	shader.use();
	shader.set_float("time", time);
	
	/*	
	glm::mat4 _model {1.0f};
	glm::mat4 _view {1.0f};
	glm::mat4 _projection {1.0f};
	glUniformMatrix4fv(glGetUniformLocation(_program, "model"), 1, GL_FALSE, glm::value_ptr(_model));
	glUniformMatrix4fv(glGetUniformLocation(_program, "view"), 1, GL_FALSE, glm::value_ptr(_view));
	glUniformMatrix4fv(glGetUniformLocation(_program, "projection"), 1, GL_FALSE, glm::value_ptr(_projection));
	*/

	glBindVertexArray(_vao);
	glDrawElements(GL_TRIANGLES, _indices.size(), GL_UNSIGNED_INT, nullptr);
	glBindVertexArray(0);
}
