#include "Entity.h"
#include <iostream>

Entity::Entity(glm::vec3 position, glm::vec3 rotation, glm::vec2 scale)
	: _position {position}
	, _rotation {rotation}
	, _scale {scale}
{

}

Entity::~Entity(void) {

}

void Entity::draw(glm::vec3 view_position, glm::mat4 projection, float delta_time) {
	std::cout << "bad draw :(" << std::endl;
}

void Entity::set_position(glm::vec3 position) {
	_position = position;
}

void Entity::set_rotation(glm::vec3 rotation) {
	_rotation = rotation;
}

glm::vec3 Entity::position() const {
	return _position;
}

glm::vec3 Entity::rotation() const {
	return _rotation;
}

void Entity::rotate_test(float delta_time) {
	_rotation.x -= 20 * delta_time;
	_rotation.y += 150 * delta_time;
	_rotation.z += 50 * delta_time;

	_position.y = sin(_rotation.z/30);
}

glm::mat4 Entity::get_model() const {
	// Matrice model pour definir la position
	// et la rotation et le scale de l'objet dans l'espace
	glm::mat4 model {1.0f};
	model = glm::translate(model, _position);
	model = glm::rotate(model, glm::radians(_rotation.x), glm::vec3(1.0f, 0.0f, 0.0f));
	model = glm::rotate(model, glm::radians(_rotation.y), glm::vec3(0.0f, 1.0f, 0.0f));
	model = glm::rotate(model, glm::radians(_rotation.z), glm::vec3(0.0f, 0.0f, 1.0f));
	model = glm::scale(model, glm::vec3{_scale, 1.0f});

	return model;
}

glm::mat4 Entity::get_view(glm::vec3 view_position) const {
	glm::mat4 view {1.0f};
	view = glm::translate(view, view_position);

	return view;
}
