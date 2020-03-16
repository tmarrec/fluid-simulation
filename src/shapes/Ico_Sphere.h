#pragma once

#include "../Shape.h"
#include "../MainWindow.h"

class Ico_Sphere : public Shape {

public:
	Ico_Sphere(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, unsigned n, MainWindow * main_window);
	~Ico_Sphere();

private:
	static std::vector<GLfloat> generate_base_vertices();
	static std::vector<GLuint> generate_base_indices();
	static glm::vec3 half_vertex(glm::vec3 v0, glm::vec3 v1);
	static void add_vertices(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, std::vector<GLfloat> &vertices);
	static void add_indices(uint64_t i, std::vector<GLuint> &indices);
	static std::tuple<std::vector<GLfloat>, std::vector<GLuint>> subdivide(std::vector<GLfloat> vertices, std::vector<GLuint> indices);
	static std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> get_geometry(unsigned n);
};
