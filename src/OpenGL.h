#pragma once

#include <GL/gl.h>
#include <GL/glext.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vector>
#include <memory>

#include "ECS.h"
#include "MainWindow.h"


class OpenGL : public QObject {
	Q_OBJECT

public:
	explicit OpenGL(unsigned int w, unsigned int h, MainWindow * main_window);
	virtual ~OpenGL(void);

	virtual void draw(void);

	unsigned short width() const;
	unsigned short height() const;
	glm::mat4 projection() const;
	void set_delta_time(float delta_time);
	float delta_time() const;
	glm::vec3 view_position() const;

	void move(uint id, char pos, float value);

public slots:
	void add_triangle();
	void add_cube();
	void add_sphere();

private:
	unsigned int _width;
	unsigned int _height;
	bool _draw_fill;

	float _delta_time;

	glm::mat4 _projection;
	glm::vec3 _view_position;

	MainWindow * _main_window;

	ECS _ecs;
};
