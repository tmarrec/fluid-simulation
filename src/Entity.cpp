#include "Entity.h"
#include <iostream>

unsigned long Entity::_next_id = 0;

Entity::Entity(std::string name, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale)
	: _position {position}
	, _rotation {rotation}
	, _scale {scale}
	, _id{_next_id}
	, _name{name + " " + std::to_string(_id)}
{
	_next_id++;	
	std::cout << "New Entity : " << _id << std::endl;
}

Entity::~Entity(void) {

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
	glm::vec3 temp = {-20*delta_time, 40*delta_time, 50*delta_time};
	set_rotation(rotation()+temp);

	auto pos = position();
	temp = {pos.x, sin(_rotation.z/35), pos.z};
	set_position(temp);
}

unsigned long Entity::id() const {
	return _id;
}

const std::string Entity::name() const {
	return _name;
}

glm::mat4 Entity::get_model() const {
	// Matrice model pour definir la position
	// et la rotation et le scale de l'objet dans l'espace
	glm::mat4 model {1.0f};
	model = glm::translate(model, _position);
	model = glm::rotate(model, glm::radians(_rotation.x), glm::vec3(1.0f, 0.0f, 0.0f));
	model = glm::rotate(model, glm::radians(_rotation.y), glm::vec3(0.0f, 1.0f, 0.0f));
	model = glm::rotate(model, glm::radians(_rotation.z), glm::vec3(0.0f, 0.0f, 1.0f));
	model = glm::scale(model, glm::vec3{_scale});

	return model;
}

glm::mat4 Entity::get_view(glm::vec3 view_position) const {
	glm::mat4 view {1.0f};
	view = glm::translate(view, view_position);
	return view;
}

