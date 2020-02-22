#include "OpenGL.h"

#include "shapes/Triangle.h"
#include "shapes/Cube.h"
#include "shapes/Sphere.h"

#include <iostream>

OpenGL::OpenGL(unsigned int w, unsigned int h, MainWindow * main_window)
		: _width{w}
		, _height{h}
		, _draw_fill{true}
		, _projection{1.0f}
		, _view_position{0.0f, 0.0f, -5.0f}
		, _main_window{main_window}
		, _ecs{}
{
	_projection = glm::perspective(glm::radians(70.0f), (float)width()/(float)height(), 0.1f, 100.0f);
	glEnable(GL_DEPTH_TEST);
	glViewport(0, 0, _width, _height);
}

OpenGL::~OpenGL(void) {

}

void OpenGL::move(uint id, char pos, float value) {
	_ecs.move(id, pos, value);
}

void OpenGL::add_triangle() {
	// TODO truc
	std::unique_ptr<Entity> s = std::unique_ptr<Triangle>(new Triangle(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f}));
	_main_window->add_item_to_QListW(s);
	_ecs.add(std::move(s));
}

void OpenGL::add_cube() {
	auto s {std::make_unique<Cube>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f})};
	//_main_window->add_item_to_QListW(s->id(), "cube");
	_ecs.add(std::move(s));
}

void OpenGL::add_sphere() {
	auto s {std::make_unique<Sphere>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f}, glm::vec2{10,10})};
	//_main_window->add_item_to_QListW(s->id(), "sphere");
	_ecs.add(std::move(s));
}

void OpenGL::draw(void) {
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if (_draw_fill) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	// Render tout ici
	_ecs.render_all(view_position(), projection(), delta_time());
	
}

unsigned short OpenGL::width() const {
	return _width;
}

unsigned short OpenGL::height() const {
	return _height;
}

glm::mat4 OpenGL::projection() const {
	return _projection;
}

float OpenGL::delta_time() const {
	return _delta_time;
}

void OpenGL::set_delta_time(float delta_time) {
	_delta_time = delta_time;
}

glm::vec3 OpenGL::view_position() const {
	return _view_position;
}
