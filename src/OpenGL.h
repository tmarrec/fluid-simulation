#pragma once

#include <GL/gl.h>
#include <GL/glext.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

#include <vector>
#include <memory>

#include "ECS.h"
#include "Camera.h"
#include "MainWindow.h"
#include "Light.h"

class OpenGL : public QObject {
	Q_OBJECT

public:
	explicit OpenGL(unsigned int w, unsigned int h, MainWindow * main_window, GLWidget * glw);
	virtual ~OpenGL(void);

	virtual void draw(void);

	unsigned short width() const;
	unsigned short height() const;
	glm::mat4 projection() const;
	void set_delta_time(float delta_time);
	float delta_time() const;
	void remove_entity(std::shared_ptr<Entity> entity);

public slots:
	void add_triangle();
	void add_cube();
	void add_sphere();
	void add_light();
	void add_model();

private:
	unsigned int _width;
	unsigned int _height;
	bool _draw_fill;

	float _delta_time;

	glm::mat4 _projection;
	std::shared_ptr<Camera> _camera;
	MainWindow * _main_window;
	ECS _ecs;
	GLWidget * _glw;
};
