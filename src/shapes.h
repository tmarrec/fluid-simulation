#pragma once

#include <vector>

#include "utils.h"

struct Shape
{
	std::vector<GLfloat> vertices;
	std::vector<GLfloat> normals;
	std::vector<GLuint> indices;
};

struct Quad
{
	std::vector<GLfloat> vertices =
	{
		-1.0f,  1.0f,  0.0f, 1.0f,
        -1.0f, -1.0f,  0.0f, 0.0f,
         1.0f, -1.0f,  1.0f, 0.0f,

        -1.0f,  1.0f,  0.0f, 1.0f,
         1.0f, -1.0f,  1.0f, 0.0f,
         1.0f,  1.0f,  1.0f, 1.0f
	};
};

struct Pyramid
{
	std::vector<GLfloat> vertices =
	{
		0.0f, 1.27201964951f, 0.0f,
		-1.0f, 0.0f, 1.0f,
		1.0f, 0.0f, 1.0f,
		-1.0f, 0.0f, -1.0f,
		1.0f, 0.0f, -1.0f,
	};
	std::vector<GLfloat> normals =
	{
		0.0f, 1.0f, 0.0f,
		-1.0f, -0.63600982475f, 1.0f,
		1.0f, -0.63600982475f, 1.0f,
		-1.0f, -0.63600982475f, -1.0f,
		1.0f, 0.63600982475f, -1.0f,
	};
	std::vector<GLuint> indices
	{
		3, 2, 1, 3, 4, 2,
		4, 0, 2,
		2, 0, 1,
		1, 0, 3,
		3, 0, 4
	};
};

struct Cube
{
	std::vector<GLfloat> vertices =
	{
		-0.5f, -0.5f, 0.5f,
		0.5f, -0.5f, 0.5f,
		0.5f, 0.5f, 0.5f,
		-0.5f, 0.5f, 0.5f,

		0.5f, 0.5f, 0.5f,
		0.5f, 0.5f, -0.5f,
		0.5f, -0.5f, -0.5f,
		0.5f, -0.5f, 0.5f,

		-0.5f, -0.5f, -0.5f,
		0.5f, -0.5f, -0.5f,
		0.5f, 0.5f, -0.5f,
		-0.5f, 0.5f, -0.5f,

		-0.5f, -0.5f, -0.5f,
		-0.5f, -0.5f, 0.5f,
		-0.5f, 0.5f, 0.5f,
		-0.5f, 0.5f, -0.5f,
						
		0.5f, 0.5f, 0.5f,
		-0.5f, 0.5f, 0.5f,
		-0.5f, 0.5f, -0.5f,
		0.5f, 0.5f, -0.5f,

		-0.5f, -0.5f, -0.5f,
		0.5f, -0.5f, -0.5f,
		0.5f, -0.5f, 0.5f,
		-0.5f, -0.5f, 0.5f,
	};
	std::vector<GLfloat> normals =
	{
		0.0f, 0.0f, 1.0f,
		0.0f, 0.0f, 1.0f,
		0.0f, 0.0f, 1.0f,
		0.0f, 0.0f, 1.0f,

		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,

		0.0f, 0.0f, -1.0f,
		0.0f, 0.0f, -1.0f,
		0.0f, 0.0f, -1.0f,
		0.0f, 0.0f, -1.0f,

		-1.0f, 0.0f, 0.0f,
		-1.0f, 0.0f, 0.0f,
		-1.0f, 0.0f, 0.0f,
		-1.0f, 0.0f, 0.0f,

		0.0f, 1.0f, 0.0f,
		0.0f, 1.0f, 0.0f,
		0.0f, 1.0f, 0.0f,
		0.0f, 1.0f, 0.0f,

		0.0f, -1.0f, 0.0f,
		0.0f, -1.0f, 0.0f,
		0.0f, -1.0f, 0.0f,
		0.0f, -1.0f, 0.0f,
	};
	std::vector<GLuint> indices
	{
		0,  1,  2,  0,  2,  3,   // Front
		4,  6,  5,  4,  7,  6,   // Right
		5,  9, 11,  11, 9,  8,    // Back
		15, 12, 13, 13, 14, 15,  // Left
		17, 16, 19, 17, 19, 18,  // Upper
		20, 21, 23, 21, 22, 23 	 // Bottom
	};
};
