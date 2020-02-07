#include "Entity.h"
#include <iostream>

Entity::Entity(glm::vec3 position, glm::vec3 rotation, glm::vec2 scale) {

}

Entity::~Entity(void) {

}

void Entity::draw(glm::vec3 view_position, glm::mat4 projection, float delta_time) {
	std::cout << "bad draw :(" << std::endl;
}
