#include "OpenGL.h"

#include "shapes/Triangle.h"
#include "shapes/BSpline_line.h"
#include "shapes/BSpline_tensor.h"
#include "shapes/Cube.h"
#include "shapes/UV_Sphere.h"
#include "shapes/Ico_Sphere.h"
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
	_camera = std::shared_ptr<Camera>(new Camera(glm::vec3{-250.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{1.0f, 1.0f, 1.0f}, 90.0f, main_window));

	_projection = glm::perspective(glm::radians(_camera->FOV()), (float)width()/(float)height(), 0.1f, 1000000.0f);
	_main_window->add_item_to_QListW(_camera);
	_ecs.add(_camera);
}

OpenGL::~OpenGL(void) {

}

// Ajoute un Triangle a la scene
void OpenGL::add_triangle() {
	_glw->make_current();
	auto s = std::shared_ptr<Triangle>(new Triangle(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{100.0f, 100.0f, 100.0f}, _main_window));
	_main_window->add_item_to_QListW(s);
	_ecs.add(s);
	_glw->done_current();
}

// Ajoute une B-Spline a la scene
void OpenGL::add_bspline() {
	_glw->make_current();
	std::vector<glm::vec3> controls = {
		{-1.5, 0.75, -1.5},
		{-0.75, -0.75, 1.5},
		{0.75, -0.75, -1.5},
		{1.5, 0.75, 1.5},
		{-1.5, 0.75, -1.5},
		{-0.75, -0.75, 1.5},
		{0.75, -0.75, -1.5},
		{1.5, 0.75, 1.5},
	};
	auto s = std::shared_ptr<BSpline_line>(new BSpline_line(4, controls, false, 0.01f, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{100.0f, 100.0f, 100.0f}, _main_window));
	_main_window->add_item_to_QListW(s);
	_ecs.add(s);
	_glw->done_current();
}

// Ajoute une B-Spline tensor a la scene
void OpenGL::add_bspline_tensor() {
	_glw->make_current();
	std::vector<std::vector<glm::vec3>> controls = {
		{
			{0, 0, 0},
			{1, 0, 0},
			{2, 0, 0}
		},
		{
			{0, 0, 1},
			{1, 3, 1},
			{2, 0, 1}
		},
		{
			{0, 0, 2},
			{1, 0, 2},
			{2, 0, 2}
		},
	};
	auto s = std::shared_ptr<BSpline_tensor>(new BSpline_tensor(3, controls, false, 0.01f, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{100.0f, 100.0f, 100.0f}, _main_window));
	_main_window->add_item_to_QListW(s);
	_ecs.add(s);
	_glw->done_current();
}

// Ajoute un Cube a la scene
void OpenGL::add_cube() {
	_glw->make_current();
	auto s = std::shared_ptr<Cube>(new Cube(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{100.0f, 100.0f, 100.0f}, _main_window));
	_main_window->add_item_to_QListW(s);
	_ecs.add(s);
	_glw->done_current();
}

// Ajoute une UV Sphere a la scene
void OpenGL::add_uv_sphere() {
	_glw->make_current();
	glm::vec2 faces {16, 16};
	auto s = std::shared_ptr<UV_Sphere>(new UV_Sphere(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{100.0f, 100.0f, 100.0f}, faces, _main_window));
	_main_window->add_item_to_QListW(s);
	_ecs.add(s);
	_glw->done_current();
}

// Ajoute une ICO Sphere a la scene
void OpenGL::add_ico_sphere() {
	_glw->make_current();
	unsigned n_subdivide = 2;
	auto s = std::shared_ptr<Ico_Sphere>(new Ico_Sphere(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{100.0f, 100.0f, 100.0f}, n_subdivide, _main_window));
	_main_window->add_item_to_QListW(s);
	_ecs.add(s);
	_glw->done_current();
}

// Ajoute une Light a la scene
void OpenGL::add_light() {
	_glw->make_current();
	float r = 1.0f;
	float g = 1.0f;
	float b = 1.0f;
	auto l = std::shared_ptr<Light>(new Light(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{2.0f, 2.0f, 2.0f}, {r, g, b}, _main_window, 1.0f));
	_main_window->add_item_to_QListW(l);
	_ecs.add(l);
	_glw->done_current();
}

// Ajoute un Model a la scene
void OpenGL::add_model(std::string model_path) {
	_glw->make_current();
	float scale = 1.0f;
	auto m = std::shared_ptr<Model>(new Model(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{scale, scale, scale}, _main_window, model_path));
	_main_window->add_item_to_QListW(m);
	_ecs.add(m);
	_glw->done_current();
}

// Change le type d'affichage opengl
void OpenGL::set_draw_fill(bool state) {
	_draw_fill = state;
}

// Clear l'ecran et draw toutes les entitées
void OpenGL::draw(void) {
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (_draw_fill) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	_ecs.render_all(_camera->view(), projection(), delta_time());
}

void OpenGL::remove_entity(std::shared_ptr<Entity> entity) {
	_ecs.remove(entity);
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

// Ajoute une lumiere qui a deja une position, une couleur et une intensitée
void OpenGL::add_light_placed(glm::vec3 position, glm::vec3 color, float intensity) {
	_glw->make_current();
	auto l = std::shared_ptr<Light>(new Light(position, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{2.0f, 2.0f, 2.0f}, color, _main_window, intensity));
	_main_window->add_item_to_QListW(l);
	_ecs.add(l);
	_glw->done_current();
}
