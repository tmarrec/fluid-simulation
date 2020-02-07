#pragma once

#include "../OpenGL.h"
#include "../Shader.h"

#include <vector>

class Triangle : public Entity {

public:
	Triangle(glm::vec3 position, glm::vec3 rotation, glm::vec2 scale);
	~Triangle();

	virtual void draw(glm::vec3 view_position, glm::mat4 projection, float delta_time);

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
	glm::vec2 _scale;

};
