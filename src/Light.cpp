#include "Light.h"
#include <iostream>

Light::Light(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
	std::string type, glm::vec3 color, MainWindow * main_window)
	: Entity(type, position, rotation, scale, main_window)
	, _color{color}
	, _type{type}
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

void Light::draw(glm::mat4 view, glm::mat4 projection, float delta_time,
				std::vector<std::shared_ptr<Entity>> lights) {
	std::ignore = view;
	std::ignore = projection;
	std::ignore = delta_time;
	std::ignore = lights;
}

void Light::set_shader(Shader shader) {
	std::cout << "Shader shouldn't be set to a light" << std::endl;
	std::ignore = shader;
}

Shader & Light::shader() {
	std::cout << "Shouldn't ask for Light shader" << std::endl;
	exit(0);
}

