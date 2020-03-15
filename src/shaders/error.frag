#version 410 core

in vec3 i_normal;
in vec3 pos_error;
out vec4 color;

uniform vec4 _color;

void main(void) {
	float err = 1 - length(pos_error);
	color = 40*vec4(err, err, err, 1.0f);
}
