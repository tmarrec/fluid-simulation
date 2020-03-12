#pragma once

#include "../Shape.h"

class Cube : public Shape {

public:
	Cube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
	~Cube();

private:

};
