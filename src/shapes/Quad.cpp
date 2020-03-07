#include "Quad.h"

#include <iostream>

Quad::Quad(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
		: Shape(
			position,
			rotation,
			scale,
			std::vector<GLfloat>{
				0.0f,		0.0f, 		0.0f, 		// Top Left
				0.5f,		-0.5f, 		0.0f, 		// Bottom Left
				-0.5f,		-0.5f, 		0.0f, 		// Bottom Right
				0.0f,		-0.5f, 		0.0f, 		// Top Right
			},
			std::vector<GLfloat>{
				0,		1,		0,
				1,		0,		1,
			},
			std::vector<GLuint>{
				0, 1, 2, 3 // Quad
			},
			"Quad",
			{1.0f, 0.5f, 1.0f},
			{
				"shaders/vert.vert",
				"shaders/frag.frag"
			}
		)
{

}

Quad::~Quad() {
	
}

