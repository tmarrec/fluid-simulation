#include "Shape.h"
#include <iostream>

Shape::Shape(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
	std::vector<GLfloat> vertices, std::vector<GLfloat> normals,
	std::vector<GLuint> indices)
	: Entity(position, rotation, scale)
	, _vertices{vertices}
	, _normals{normals}
	, _indices{indices}
{
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
			vertices.size()*sizeof(GLfloat),
			vertices.data(),
			GL_STATIC_DRAW
		);

		// 4. Then set our vertex attributes pointers
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(0);

		// 5. Copy our normals array in a buffer for OpenGL to use
		glBindBuffer(GL_ARRAY_BUFFER, _nbo);
		glBufferData(
			GL_ARRAY_BUFFER,
			normals.size()*sizeof(GLfloat),
			normals.data(),
			GL_STATIC_DRAW
		);

		// 6. Copy our vertices array in a buffer for OpenGL to use
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(1);

		// 7. Copy our index array in a element buffer for OpenGL to use
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _ebo);
		glBufferData(
			GL_ELEMENT_ARRAY_BUFFER,
			indices.size()*sizeof(GLfloat),
			indices.data(),
			GL_STATIC_DRAW
		);

	// 8. Unbind the VAO
	glBindVertexArray(0);

}

Shape::~Shape(void) {

}

std::vector<GLfloat> Shape::vertices() const {
	return _vertices;
}

std::vector<GLfloat> Shape::normals() const {
	return _normals;
}

std::vector<GLuint> Shape::indices() const {
	return _indices;
}

void Shape::draw_vertex() const {
	glBindVertexArray(_vao);
	std::cout << indices().size() << std::endl;
	glDrawElements(GL_TRIANGLES, indices().size(), GL_UNSIGNED_INT, nullptr);
	glBindVertexArray(0);
}
