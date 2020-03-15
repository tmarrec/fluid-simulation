#include "Camera.h"
#include <iostream>

Camera::Camera(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
	float FOV, MainWindow * main_window)
	: Entity(position, rotation, scale, main_window)
	, _FOV{FOV}
	, _front{glm::vec3{0.0f, 0.0f, 1.0f}}
	, _up{glm::vec3{0.0f, 1.0f, 0.0f}}
	, _yaw{0}
	, _pitch{0}
{

}

Camera::~Camera(void) {

}

Entity_Type Camera::type() {
	return CAMERA;
}

float Camera::FOV() const {
	return _FOV;
}

float Camera::yaw() const {
	return _yaw;
}

float Camera::pitch() const {
	return _pitch;
}

glm::vec3 Camera::front() const {
	return _front;
}

glm::vec3 Camera::up() const {
	return _up;
}

void Camera::set_yaw(float yaw) {
	_yaw = yaw;
}

void Camera::set_pitch(float pitch) {
	_pitch = pitch;
}

glm::mat4 Camera::view() const {
	return glm::lookAt(position(), position()+_front, _up);
}

void Camera::set_front(glm::vec3 front) {
	_front = front;
}

void Camera::draw(glm::mat4 view, glm::mat4 projection, float delta_time,
				std::vector<std::shared_ptr<Entity>> lights) {
	std::ignore = view;
	std::ignore = projection;
	std::ignore = delta_time;
	std::ignore = lights;
}

void Camera::set_shader(Shader shader) {
	std::cout << "Shader shouldn't be set to a camera" << std::endl;
	std::ignore = shader;
}

Shader & Camera::shader() {
	std::cout << "Shouldn't ask for Camera shader" << std::endl;
	exit(0);
}

