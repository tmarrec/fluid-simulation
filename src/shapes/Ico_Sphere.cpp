#include "Ico_Sphere.h"

#include <iostream>

Ico_Sphere::Ico_Sphere(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, glm::vec2 faces, MainWindow * main_window)
		: Shape(
			position,
			rotation,
			scale,
			get_geometry(faces),
			{RND(), RND(), RND()},
			{
				"shaders/vert.vert",
				"shaders/frag.frag"
			},
			main_window
		)
{
	set_geometry();
}

std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> Ico_Sphere::get_geometry(glm::vec2 faces) {
	auto vertices = generate_vertices(faces);
	return {vertices, vertices, generate_indices(faces)};
}

std::vector<GLfloat> Ico_Sphere::subdivide() {

}

std::vector<GLfloat> Ico_Sphere::generate_vertices(glm::vec2 faces) {
	std::vector<GLfloat> vertices(12*3);
	float a = 0.52573111211;
	float b = 0.85065080835;

	vertices[0] = a; vertices[1] = b; vertices[2] = 0;
	vertices[3] = b; vertices[4] = 0; vertices[5] = a;
	vertices[6] = 0; vertices[7] = a; vertices[8] = b;

	vertices[9] = -a; vertices[10] = b; vertices[11] = 0;
	vertices[12] = 0; vertices[13] = a; vertices[14] = -b;
	vertices[15] = b; vertices[16] = 0; vertices[17] = -a;

	vertices[18] = a; vertices[19] = -b; vertices[20] = 0;
	vertices[21] = 0; vertices[22] = -a; vertices[23] = b;
	vertices[24] = -b; vertices[25] = 0; vertices[26] = a;

	vertices[27] = -b; vertices[28] = 0; vertices[29] = -a;
	vertices[30] = 0; vertices[31] = -a; vertices[32] = -b;
	vertices[33] = -a; vertices[34] = -b; vertices[35] = 0;

	return vertices;
}

std::vector<GLuint> Ico_Sphere::generate_indices(glm::vec2 faces) {
	std::vector<GLuint> indices(20*3);
	indices[0] = 0; 	indices[1] = 2; 	indices[2] = 1;
	indices[3] = 0; 	indices[4] = 1; 	indices[5] = 5;
	indices[6] = 0; 	indices[7] = 5; 	indices[8] = 4;
	indices[9] = 0; 	indices[10] = 4; 	indices[11] = 3;
	indices[12] = 0; 	indices[13] = 3; 	indices[14] = 2;

	indices[15] = 11; 	indices[16] = 10; 	indices[17] = 6;
	indices[18] = 11; 	indices[19] = 6; 	indices[20] = 7;
	indices[21] = 11; 	indices[22] = 7; 	indices[23] = 8;
	indices[24] = 11; 	indices[25] = 8; 	indices[26] = 9;
	indices[27] = 11; 	indices[28] = 9; 	indices[29] = 10;

	indices[30] = 1; 	indices[31] = 7; 	indices[32] = 6;
	indices[33] = 5; 	indices[34] = 6; 	indices[35] = 10;
	indices[36] = 4; 	indices[37] = 10; 	indices[38] = 9;
	indices[39] = 3; 	indices[40] = 9; 	indices[41] = 8;
	indices[42] = 2; 	indices[43] = 8; 	indices[44] = 7;

	indices[45] = 10; 	indices[46] = 4; 	indices[47] = 5;
	indices[48] = 9; 	indices[49] = 3; 	indices[50] = 4;
	indices[51] = 8; 	indices[52] = 2; 	indices[53] = 3;
	indices[54] = 7; 	indices[55] = 1; 	indices[56] = 2;
	indices[57] = 6; 	indices[58] = 5; 	indices[59] = 1;
	return indices;
}

Ico_Sphere::~Ico_Sphere() {
	
}
