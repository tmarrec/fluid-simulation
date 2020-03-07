#pragma once

#include "../Shape.h"

class Cube : public Shape {

public:
	Cube(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
	~Cube();

	void draw(glm::vec3 view_position, glm::mat4 projection, float delta_time, std::vector<std::shared_ptr<Entity>> lights) override;


private:
	void apply_color();
};
