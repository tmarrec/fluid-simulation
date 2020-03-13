#include "Light.h"
#include <iostream>

Light::Light(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
	std::string type, glm::vec3 color, MainWindow * main_window)
	: Entity(type, position, rotation, scale, main_window)
	, _type{type}
	, _color{color}
{

}

Light::~Light(void) {

}

const std::string Light::type() const {
	return _type;
}

const glm::vec3 Light::color() const {
	return _color;
}

void Light::draw(glm::mat4 view_position, glm::mat4 projection, float delta_time,
				std::vector<std::shared_ptr<Entity>> lights) {

}

