#include "ECS.h"
#include <iostream>

ECS::ECS() {

}

ECS::~ECS(void) {

}

void ECS::add(std::unique_ptr<Entity> entity) {
	_entities.push_back(std::move(entity));
}

void ECS::move(uint id, char pos, float value) {
	for (auto & e : _entities) {
		if (e->id() == id) {
			auto p = e->position();
			if (pos == 'x') {
				e->set_position(glm::vec3{0.05f*value, p.y, p.z});
			} else if (pos == 'y') {
				e->set_position(glm::vec3{p.x, 0.05f*value, p.z});
			} else if (pos == 'z') {
				e->set_position(glm::vec3{p.x, p.y, 0.05f*value});
			}	
		}
	}
}

void ECS::render_all(glm::vec3 view_position, glm::mat4 projection, float delta_time) {
	for (const auto & e : _entities) {
		e->draw(view_position, projection, delta_time);
	}
}

