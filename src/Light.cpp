#include "Light.h"
#include <iostream>

Light::Light(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
	glm::vec3 color, MainWindow * main_window, float intensity)
	: Entity(position, rotation, scale, main_window)
	, _color{color}
	, _intensity{intensity}
{

}

Light::~Light(void) {

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

Entity_Type Light::type() {
	return LIGHT;
}

float Light::intensity() const {
	return _intensity;
}
