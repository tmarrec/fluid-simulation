#include "ECS.h"
#include <iostream>

ECS::ECS() {

}

ECS::~ECS(void) {

}

void ECS::add(std::shared_ptr<Entity> entity) {
	_entities.push_back(std::move(entity));
}

void ECS::render_all(glm::vec3 view_position, glm::mat4 projection, float delta_time) {
	for (const auto & e : _entities) {
		e->draw(view_position, projection, delta_time);
	}
}

