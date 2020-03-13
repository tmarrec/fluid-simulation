#pragma once

#include "../Shape.h"
#include "../MainWindow.h"

class Triangle : public Shape {

public:
	Triangle(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale, MainWindow * main_windown);
	~Triangle();

private:

};
