#pragma once

#include "../OpenGL.h"

#include "../Shader.h"

#include <vector>

class Triangle : public OpenGL {

public:
	explicit Triangle(int w, int h);
	~Triangle() override;

	void draw() override;
	float test = 0;


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

};
