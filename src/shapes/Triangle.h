#pragma once

#include "../OpenGL.h"

#include "../Shader.h"

#include <vector>

class Triangle : public OpenGL {

public:
	explicit Triangle(int w, int h);
	~Triangle() override;

	void draw() override;

	void set_position(glm::vec3 position);
	void set_rotation(glm::vec3 rotation);
	glm::vec3 position() const;
	glm::vec3 rotation() const;

private:
	std::vector<GLfloat> 	_vertices;
	std::vector<GLfloat> 	_normals;
	std::vector<GLuint>		_indices;

	// Vertex Array Buffer
	GLuint _vao;
	// Vertex Buffer Object
	GLuint _vbo;
	// Normal Buffer
	GLuint _nbo;
	// Face Buffer
	GLuint _ebo;

	// Shader Program
	GLuint _program;

	Shader _shader;

	glm::vec4 _color;
	void apply_color();

	glm::vec3 _position;
	glm::vec3 _rotation;

};
