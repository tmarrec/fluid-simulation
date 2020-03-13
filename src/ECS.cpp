#include "ECS.h"
#include <iostream>
#include <algorithm>

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

void ECS::remove(std::shared_ptr<Entity> entity) {
	// cf le TODO au dessus
	if (entity->name().find("Light") != std::string::npos) {
		_lights.erase(std::remove(_lights.begin(), _lights.end(), entity), _lights.end());
	}
	_entities.erase(std::remove(_entities.begin(), _entities.end(), entity), _entities.end());
}


std::vector<std::shared_ptr<Entity>> ECS::lights() const {
	return _lights;
}

void ECS::render_all(glm::mat4 view_position, glm::mat4 projection, float delta_time) {
	for (const auto & e : _entities) {
		e->draw(view_position, projection, delta_time, lights());
	}
}

