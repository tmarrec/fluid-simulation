#include "Camera.h"
#include <iostream>

Camera::Camera(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
	std::string type, float FOV, MainWindow * main_window)
	: Entity(type, position, rotation, scale, main_window)
	, _type{type}
	, _FOV{FOV}
{

}

Camera::~Camera(void) {

}

const std::string Camera::type() const {
	return _type;
}

float Camera::FOV() const {
	return _FOV;
}

void Camera::draw(glm::vec3 view_position, glm::mat4 projection, float delta_time,
				std::vector<std::shared_ptr<Entity>> lights) {

}

