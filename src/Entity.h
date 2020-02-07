#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class Entity {

public:
	Entity(glm::vec3 position, glm::vec3 rotation, glm::vec2 scale);
	~Entity(void);

	virtual void draw(glm::vec3 view_position, glm::mat4 projection, float delta_time);

private:

};
