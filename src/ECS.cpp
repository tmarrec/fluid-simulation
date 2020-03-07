#include "ECS.h"
#include <iostream>

ECS::ECS() {

}

ECS::~ECS(void) {

}

void ECS::add(std::shared_ptr<Entity> entity) {
	// TODO très moche, devrais plutot faire fonction virtuelle "type()" sur 
	// Entity.h pour recuperer le type de chaque entitées
	if (entity->name().find("Light") != std::string::npos) {
		_lights.push_back(entity);
	}
	_entities.push_back(entity);
}


std::vector<std::shared_ptr<Entity>> ECS::lights() const {
	return _lights;
}

void ECS::render_all(glm::vec3 view_position, glm::mat4 projection, float delta_time) {
	for (const auto & e : _entities) {
		e->draw(view_position, projection, delta_time, lights());
	}
}

