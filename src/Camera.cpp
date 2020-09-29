#include "Camera.h"
#include <iostream>

Camera::Camera(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
	float FOV, MainWindow * main_window)
	: Entity(position, rotation, scale, main_window)
	, _FOV{FOV}
	, _front{glm::vec3{1.0f, 0.0f, 0.0f}}
	, _up{glm::vec3{0.0f, 1.0f, 0.0f}}
	, _yaw{0}
	, _pitch{0}
	, _speed{8}
	, _move_front{false}
	, _move_back{false}
	, _move_left{false}
	, _move_right{false}
	, _move_up{false}
	, _move_down{false}
{

}

Camera::~Camera(void) {

}

Entity_Type Camera::type() {
	return CAMERA;
}

float Camera::speed() const {
	return _speed;
}

void Camera::set_speed(float speed) {
	_speed = speed;
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

void Camera::set_move_front(bool state) {
	_move_front = state;
}

void Camera::set_move_back(bool state) {
	_move_back = state;
}

void Camera::set_move_left(bool state) {
	_move_left = state;
}

void Camera::set_move_right(bool state) {
	_move_right = state;
}

void Camera::set_move_up(bool state) {
	_move_up = state;
}

void Camera::set_move_down(bool state) {
	_move_down = state;
}

bool Camera::move_front() const {
	return _move_front;
}

bool Camera::move_back() const {
	return _move_back;
}

bool Camera::move_left() const {
	return _move_left;
}

bool Camera::move_right() const {
	return _move_right;
}

bool Camera::move_up() const {
	return _move_up;
}

bool Camera::move_down() const {
	return _move_down;
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

