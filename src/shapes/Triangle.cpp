#include "Triangle.h"

#include <iostream>

Triangle::Triangle(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
		: Shape(
			position,
			rotation,
			scale,
			std::vector<GLfloat>{
				0.0f,		0.5f, 		0.0f, 		// Top Center
				0.5f,		-0.5f, 		0.0f, 		// Bottom Left
				-0.5f,		-0.5f, 		0.0f, 		// Bottom Right
			},
			std::vector<GLfloat>{
				1,		1,		1,
				1,		1,		1,
				1,		1,		1,
			},
			std::vector<GLuint>{
				0, 2, 1, // Triangle
			},
			"Triangle",
			{1.0f, 0.5f, 1.0f},
			{
				"shaders/vert.vert",
				"shaders/frag.frag"
			}
		)
{

}

Triangle::~Triangle() {
	
}

