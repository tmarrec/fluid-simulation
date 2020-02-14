#include "ECS.h"
#include <iostream>

ECS::ECS() {

}

ECS::~ECS(void) {

}

void ECS::add(std::unique_ptr<Entity> entity) {
	_entities.push_back(std::move(entity));
}

void ECS::get(uint id) {
	for (auto & e : _entities) {
		if (e->id() == id) {
			e->set_position(glm::vec3{0.0f, 0.0f, 2.0f});
		}
	}
}

void ECS::render_all(glm::vec3 view_position, glm::mat4 projection, float delta_time) {
	for (const auto & e : _entities) {
		e->draw(view_position, projection, delta_time);
	}
}

