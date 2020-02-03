#include "OpenGL.h"

#include <iostream>

OpenGL::OpenGL(int w, int h)
		: _width{w}
		, _height{h}
		, _draw_fill{true} {
	glEnable(GL_DEPTH_TEST);
	glViewport(0, 0, _width, _height);
}

OpenGL::~OpenGL(void) {

}

void OpenGL::draw(void) {
	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if (_draw_fill) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
}
