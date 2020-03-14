#include "Entity.h"
#include <iostream>

#include "MainWindow.h"

unsigned long Entity::_next_id = 0;

Entity::Entity(std::string name, glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
		MainWindow * main_window)
	: _position {position}
	, _rotation {rotation}
	, _scale {scale}
	, _id{_next_id}
	, _name{name + " " + std::to_string(_id)}
	, _main_window{main_window}
{
	_next_id++;	
}

Entity::~Entity(void) {

}

void Entity::set_position(glm::vec3 position) {
	_position = position;
	_main_window->update_slide_position(_position, _id);
}

void Entity::set_rotation(glm::vec3 rotation) {
	_rotation = rotation;

	_rotation.x = std::fmod(rotation.x, 360.0f);
	if (_rotation.x < 0.0f) {
		_rotation.x += 360.0f;
	}
	
	_rotation.y = std::fmod(rotation.y, 360.0f);
	if (_rotation.y < 0.0f) {
		_rotation.y += 360.0f;
	}

	_rotation.z = std::fmod(rotation.z, 360.0f);
	if (_rotation.z < 0.0f) {
		_rotation.z += 360.0f;
	}

	_main_window->update_slide_rotation(_rotation, _id);
}

void Entity::set_scale(glm::vec3 scale) {
	_scale = scale;
	_main_window->update_slide_scale(_scale, _id);
}

glm::vec3 Entity::position() const {
	return _position;
}

glm::vec3 Entity::rotation() const {
	return _rotation;
}

glm::vec3 Entity::scale() const {
	return _scale;
}

void Entity::rotate_test(float delta_time) {
	glm::vec3 temp = {-30*delta_time, 50*delta_time, 70*delta_time};
	set_rotation(rotation()+temp);
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

