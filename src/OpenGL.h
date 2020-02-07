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

	unsigned short width() const;
	unsigned short height() const;
	glm::mat4 projection() const;
	void set_delta_time(float delta_time);
	float delta_time() const;
	glm::vec3 view_position() const;

private:
	unsigned short _width;
	unsigned short _height;
	bool _draw_fill;

	float _delta_time;

	glm::mat4 _projection;
	glm::vec3 _view_position;
};
