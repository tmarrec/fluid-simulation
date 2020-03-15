#include "ECS.h"
#include <iostream>
#include <algorithm>

ECS::ECS() {

}

ECS::~ECS(void) {

}

// Ajoute une entité aux différents vecteurs d'entitées
void ECS::add(std::shared_ptr<Entity> entity) {
	if (entity->type() == LIGHT) {
		_lights.push_back(entity);
	}
	_entities.push_back(entity);
}

// Enleve une entité aux différents vecteurs auquel elle appartient
void ECS::remove(std::shared_ptr<Entity> entity) {
	if (entity->type() == LIGHT) {
		_lights.erase(std::remove(_lights.begin(), _lights.end(), entity), _lights.end());
	}
	_entities.erase(std::remove(_entities.begin(), _entities.end(), entity), _entities.end());
}


std::vector<std::shared_ptr<Entity>> ECS::lights() const {
	return _lights;
}

// Affiche toutes les entitées de _entities
void ECS::render_all(glm::mat4 view_position, glm::mat4 projection, float delta_time) {
	for (const auto & e : _entities) {
		e->draw(view_position, projection, delta_time, lights());
	}
}

