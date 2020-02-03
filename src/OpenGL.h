#pragma once

#include <GL/gl.h>
#include <GL/glext.h>

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
