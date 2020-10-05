#pragma once

// TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO
// Create multiple messages class, this one
// is trash af


#include <vector>
#include <iostream>
#include <memory>
#include "Shader.h"

#ifdef __APPLE__
	#include <OpenGL/gl.h>
#else
	#include <GL/gl.h>
#endif


class System;

enum Type
{
	 TEST
	,HELLO
	,HELLO_ACK
	,INIT_GL
	,DRAW
	,ASK_CAMERA_INFOS_FOR_DRAW
	,ASK_ENTITIES_DRAW
	,RESIZE_GL
	,INIT_DRAWABLE
	,FREE_DRAWABLE
};

class Message
{
public:

	Message(Type type)
		: _type{type}
	{};

	Message(Type type, System * system)
		: _type{type}
		, _system{system}
	{};

	Message(Type type, int width, int height)
		: _type{type}
		, _width{width}
		, _height{height}
	{};

	Message(Type type, std::vector<GLfloat> vertices,
		std::vector<GLfloat> normals,
		std::vector<GLuint>  indices,
		GLuint VAO,
		GLuint VBO,
		GLuint NBO,
		GLuint EBO,
		glm::vec3 color,
		std::shared_ptr<Shader> shader)
		: _type{type}
		, _vertices{vertices}
		, _normals{normals}
		, _indices{indices}
		, _VAO{VAO}
		, _VBO{VBO}
		, _NBO{NBO}
		, _EBO{EBO}
		, _color{color}
		, _shader{shader}
	{};

	Message(Type type, std::vector<GLfloat> vertices,
		std::vector<GLfloat> normals,
		std::vector<GLuint>  indices,
		GLuint VAO,
		GLuint VBO,
		GLuint NBO,
		GLuint EBO,
		glm::vec3 color,
		std::shared_ptr<Shader> shader,
		glm::vec3 position,
		glm::vec3 rotation,
		glm::vec3 scale
		)
		: _type{type}
		, _vertices{vertices}
		, _normals{normals}
		, _indices{indices}
		, _VAO{VAO}
		, _VBO{VBO}
		, _NBO{NBO}
		, _EBO{EBO}
		, _color{color}
		, _shader{shader}
		, _position{position}
		, _rotation{rotation}
		, _scale{scale}
	{};

	const Type & type() const;
	System * system() const;
	int width() const;
	int height() const;


	// Ugly but will refactor later anyway 
	const Type _type;
	System * _system = nullptr;
	int _width = 0;
	int _height = 0;

	std::vector<GLfloat> _vertices;
	std::vector<GLfloat> _normals;
	std::vector<GLuint>  _indices;

	GLuint _VAO;
	GLuint _VBO;
	GLuint _NBO;
	GLuint _EBO;

	glm::vec3 _color;
	std::shared_ptr<Shader> _shader;

	glm::vec3 _position;
	glm::vec3 _rotation;
	glm::vec3 _scale;


private:
};
