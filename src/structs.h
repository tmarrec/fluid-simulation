#pragma once

struct Shape
{
	std::vector<GLfloat> vertices;
	std::vector<GLfloat> normals;
	std::vector<GLuint> indices;
};

struct GLObjects
{
	GLuint VAO;
	GLuint VBO;
	GLuint NBO;
	GLuint EBO;
};
