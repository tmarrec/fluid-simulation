#version 450 core

in vec3 normal;
in vec3 pos;
out vec4 color;

uniform vec4 _color;

void main(void) {
	float err = 1 - length(pos);
	color = 40*vec4(err, err, err, 1.0f);
}
