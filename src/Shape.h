#pragma once

#include "OpenGL.h"
#include "Shader.h"

class Shape : public Entity {

public:
	Shape(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
		std::vector<GLfloat> vertices, std::vector<GLfloat> normals,
		std::vector<GLuint> indices, std::string type);
	~Shape(void);
	
	std::vector<GLuint> indices() const;

	void draw_vertex();
	void draw_vertex_quads();
	
	const std::string type() const;

private:
	std::vector<GLfloat> 	_vertices;
	std::vector<GLfloat> 	_normals;
	std::vector<GLuint>		_indices;
	
	// Vertex Array Buffer
	GLuint _VAO;
	
	const std::string _type;
};
