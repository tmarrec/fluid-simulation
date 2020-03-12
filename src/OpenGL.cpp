#include "OpenGL.h"

#include "shapes/Triangle.h"
#include "shapes/Cube.h"
#include "shapes/Sphere.h"
#include "shapes/Model.h"

#include <iostream>

OpenGL::OpenGL(unsigned int w, unsigned int h, MainWindow * main_window, GLWidget * glw)
		: _width{w}
		, _height{h}
		, _draw_fill{true}
		, _projection{1.0f}
		, _main_window{main_window}
		, _ecs{}
		, _glw{glw}
{

	//auto c = std::shared_ptr<Camera>(&_camera);
	_camera = std::shared_ptr<Camera>(new Camera(glm::vec3{0.0f, -50.0f, -100.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f}, "Camera", 70.0f, main_window));

	_projection = glm::perspective(glm::radians(_camera->FOV()), (float)width()/(float)height(), 0.1f, 1000.0f);
	//_main_window->add_item_to_QListW(_camera);
	//_ecs.add(_camera);
}

OpenGL::~OpenGL(void) {

}

void OpenGL::add_triangle() {
	_glw->make_current();
	auto s = std::shared_ptr<Triangle>(new Triangle(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f}));
	_main_window->add_item_to_QListW(s);
	_ecs.add(s);
	_glw->done_current();
}

void OpenGL::add_cube() {
	_glw->make_current();
	auto s = std::shared_ptr<Cube>(new Cube(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f}, _main_window));
	_main_window->add_item_to_QListW(s);
	_ecs.add(s);
	_glw->done_current();
}

void OpenGL::add_sphere() {
	_glw->make_current();
	auto s = std::shared_ptr<Sphere>(new Sphere(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{2.0f, 2.0f, 2.0f}, glm::vec2{16, 16}));
	_main_window->add_item_to_QListW(s);
	_ecs.add(s);
	_glw->done_current();
}

void OpenGL::add_light() {
	_glw->make_current();
	float r = 1.0f;
	float g = 1.0f;
	float b = 1.0f;
	auto l = std::shared_ptr<Light>(new Light(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{2.0f, 2.0f, 2.0f}, "Light", {r, g, b}, _main_window));
	_main_window->add_item_to_QListW(l);
	_ecs.add(l);
	_glw->done_current();
}

void OpenGL::add_model() {
	_glw->make_current();
	float scale = 0.7f;
	auto m = std::shared_ptr<Model>(new Model(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{scale, scale, scale}));
	_main_window->add_item_to_QListW(m);
	//m->test();
	_ecs.add(m);
	_glw->done_current();
}

void OpenGL::draw(void) {
	glClearColor(0.0f, 0.0f, 0.2f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (_draw_fill) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}

	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	_ecs.render_all(_camera->position(), projection(), delta_time());
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

