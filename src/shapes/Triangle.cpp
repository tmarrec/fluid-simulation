#include "Triangle.h"
#include "../OpenGL.h"
#include "../Shader.h"

#include <iostream>
#include <string>
#include <fstream>
#include <streambuf>

Triangle::Triangle(int w, int h)
		: OpenGL(w, h)
{
	std::cout << "Triangle construct" << std::endl;
	_vertices = {
		0.5f,		0.5f, 		0.0f, 		// Top Right
		0.5f,		-0.5f, 		0.0f, 		// Bottom Right
		-0.5f,		-0.5f, 		0.0f, 		// Bottom Left
		-0.5f,		0.5f, 		0.0f, 		// Bottom Right
	};

	_normals = {
		0,		1,		0,
		0,		0,		0,
		0,		0,		0,
	};

	_indices = {
		0, 1, 3, // First Triangle
		1, 2, 3  // Second Triangle
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
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)0);
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
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)0);
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

	Shader test {
			"../src/shaders/triangle.vert",
			"../src/shaders/triangle.frag"
	};
	test.use();


}

Triangle::~Triangle() {
	
}

void Triangle::draw() {
	std::cout << "draw triangle" << std::endl;
	OpenGL::draw();

	glUseProgram(_program);
	glBindVertexArray(_vao);
	glDrawElements(GL_TRIANGLES, _indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}
