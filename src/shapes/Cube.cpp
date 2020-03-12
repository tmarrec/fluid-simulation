#include "Cube.h"


#include <iostream>

Cube::Cube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_window)
		: Shape(
			position,
			rotation,
			scale,
			{   // Vertices
				std::vector<GLfloat>{
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
				},
				// Normals
				std::vector<GLfloat>{
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
				},
				// Indices
				std::vector<GLuint>{
					0,  1,  2,  0,  2,  3,   //front
					4,  5,  6,  4,  6,  7,   //right
					8,  9,  10, 8,  10, 11,  //back
					12, 13, 14, 12, 14, 15,  //left
					16, 17, 18, 16, 18, 19,  //upper
					20, 21, 22, 20, 22, 23 	 //bottom
				}
			},
			"Cube",
			{1.0f, 1.0f, 0.0f},
			{
				"shaders/vert.vert",
				"shaders/frag.frag"
			},
			main_window
		)
{
	set_geometry();
}

Cube::~Cube() {
	
}
