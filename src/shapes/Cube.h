#pragma once

#include "../Shape.h"
#include "../MainWindow.h"

class Cube : public Shape {

public:
	Cube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_window);
	~Cube();

private:

};
