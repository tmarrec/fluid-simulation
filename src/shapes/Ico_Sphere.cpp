#include "Ico_Sphere.h"

#include <iostream>

Ico_Sphere::Ico_Sphere(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, unsigned n, MainWindow * main_window)
		: Shape(
			position,
			rotation,
			scale,
			get_geometry(n),
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

Ico_Sphere::~Ico_Sphere() {
	
}

std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> Ico_Sphere::get_geometry(unsigned n) {
	std::tuple<std::vector<GLfloat>, std::vector<GLuint>> subdivided = {generate_base_vertices(), generate_base_indices()};
	// Applique n subdivisions
	for (unsigned i = 0; i < n; ++i) {
		subdivided = subdivide(std::get<0>(subdivided), std::get<1>(subdivided));
	}
	return {std::get<0>(subdivided), std::get<0>(subdivided), std::get<1>(subdivided)};
}

// Récupere le vertex sur la moitié du sommet entre v0 et v1
glm::vec3 Ico_Sphere::half_vertex(glm::vec3 v0, glm::vec3 v1) {
	glm::vec3 new_v = {v0[0]+v1[0], v0[1]+v1[1], v0[2]+v1[2]};
	float scale = 1.0f / sqrtf(new_v[0]*new_v[0] + new_v[1]*new_v[1] + new_v[2]*new_v[2]);
	new_v *= scale;
	return new_v;
}

// Ajoute le triangle aux vertices de la sphere
void Ico_Sphere::add_vertices(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, std::vector<GLfloat> &vertices) {
	vertices.push_back(v0[0]);
	vertices.push_back(v0[1]);
	vertices.push_back(v0[2]);
	vertices.push_back(v1[0]);
	vertices.push_back(v1[1]);
	vertices.push_back(v1[2]);
	vertices.push_back(v2[0]);
	vertices.push_back(v2[1]);
	vertices.push_back(v2[2]);
}

void Ico_Sphere::add_indices(uint64_t i, std::vector<GLuint> &indices) {
	indices.push_back(i);
	indices.push_back(i+1);
	indices.push_back(i+2);
}

std::tuple<std::vector<GLfloat>, std::vector<GLuint>> Ico_Sphere::subdivide(std::vector<GLfloat> vertices, std::vector<GLuint> indices) {
	std::vector<GLfloat> temp_vertices = vertices;
	std::vector<GLuint> temp_indices = indices;
	vertices.clear();
	indices.clear();
	glm::vec3 new_v0, new_v1, new_v2;
	uint64_t index = 0;
	
	for (uint64_t i = 0; i < temp_indices.size(); i += 3) {
		// Récupere les 3 points du triangle
		glm::vec3 v0 = glm::vec3(temp_vertices[temp_indices[i]*3], temp_vertices[temp_indices[i]*3+1], temp_vertices[temp_indices[i]*3+2]);
		glm::vec3 v1 = glm::vec3(temp_vertices[temp_indices[i+1]*3], temp_vertices[temp_indices[i+1]*3+1], temp_vertices[temp_indices[i+1]*3+2]);
		glm::vec3 v2 = glm::vec3(temp_vertices[temp_indices[i+2]*3], temp_vertices[temp_indices[i+2]*3+1], temp_vertices[temp_indices[i+2]*3+2]);

		// Prend la moitié de chaque coté du triangle trouvé
		new_v0 = half_vertex(v0, v1);
		new_v1 = half_vertex(v1, v2);
		new_v2 = half_vertex(v0, v2);

		/* Ajoute les 4 nouveaux triangles
		//        v0
		//       /  \         <- A
		//      /    \
		//   new_v0-new_v2
		//    / \    /  \     <- B, C, D respectivement
		//   /   \  /    \
		//  v1---new_v1---v2
		*/
		add_vertices(v0, new_v0, new_v2, vertices); // A
		add_vertices(new_v0, v1, new_v1, vertices); // B
		add_vertices(new_v0, new_v1, new_v2, vertices); // C
		add_vertices(new_v2, new_v1, v2, vertices); // D

		add_indices(index, indices);
		add_indices(index+3, indices);
		add_indices(index+6, indices);
		add_indices(index+9, indices);
		index += 12;
	}
	return {vertices, indices};
}

// Vertex de l'icosaedre
std::vector<GLfloat> Ico_Sphere::generate_base_vertices() {
	std::vector<GLfloat> vertices(12*3);
	float a = 0.52573111211f;
	float b = 0.85065080835f;

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

// Indices de l'icosaedre
std::vector<GLuint> Ico_Sphere::generate_base_indices() {
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

