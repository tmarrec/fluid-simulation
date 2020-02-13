#include "OpenGL.h"

#include "shapes/Triangle.h"
#include "shapes/Cube.h"
#include "shapes/Sphere.h"

#include <iostream>

OpenGL::OpenGL(int w, int h, MainWindow * main_window)
		: _width{w}
		, _height{h}
		, _draw_fill{true}
		, _projection{1.0f}
		, _view_position{0.0f, 0.0f, -5.0f}
		, _ecs{}
		, _main_window{main_window}
{

	_projection = glm::perspective(glm::radians(70.0f), (float)width()/(float)height(), 0.1f, 100.0f);
	glEnable(GL_DEPTH_TEST);
	glViewport(0, 0, _width, _height);


	// Tests
	auto t1 {std::make_unique<Triangle>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f})};
	auto t2 {std::make_unique<Cube>(glm::vec3{3.0f, 0.0f, 0.0f}, glm::vec3{5.0f, 80.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f})};
	auto t3 {std::make_unique<Sphere>(glm::vec3{-3.0f, 0.0f, 0.0f}, glm::vec3{5.0f, 80.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f}, glm::vec2{5, 5})};

	_ecs.add(std::move(t1));
	_ecs.add(std::move(t2));
	_ecs.add(std::move(t3));
	// Fin Test
}

OpenGL::~OpenGL(void) {

}

void OpenGL::add_shape(std::string shape) {
	// TODO a changer avec une enum plutot que string
	if (shape == "triangle") {
		auto s {std::make_unique<Triangle>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f})};
		std::cout << s->id() << std::endl;
		_main_window->add_item_to_QListW(shape+std::to_string(s->id()));
		_ecs.add(std::move(s));
	} else if (shape == "cube") {
		auto s {std::make_unique<Cube>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f})};
		_main_window->add_item_to_QListW(shape+std::to_string(s->id()));
		_ecs.add(std::move(s));
	} else {
		auto s {std::make_unique<Sphere>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f}, glm::vec2{5, 5})};
		_main_window->add_item_to_QListW(shape+std::to_string(s->id()));
		_ecs.add(std::move(s));
	}

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
