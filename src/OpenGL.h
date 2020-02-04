#pragma once

#include <GL/gl.h>
#include <GL/glext.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class OpenGL {

public:
	explicit OpenGL(int w, int h);
	virtual ~OpenGL(void);

	virtual void draw(void);

private:
	unsigned short _width;
	unsigned short _height;
	bool _draw_fill;

};
