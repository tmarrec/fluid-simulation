#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class Entity {

public:
	Entity(void);
	~Entity(void);

	virtual void draw(glm::vec3 view_position, glm::mat4 projection, float delta_time);

private:

};
