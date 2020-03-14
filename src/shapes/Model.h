#pragma once

#include "../Shape.h"
#include "../MainWindow.h"

class Model : public Shape {

public:
	Model(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_window,
			std::string model_path);
	~Model();

private:
	static std::tuple<std::vector<GLfloat>, std::vector<GLfloat>, std::vector<GLuint>> get_geometry(
			std::string model_path);
};
