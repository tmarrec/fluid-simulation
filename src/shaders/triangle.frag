#version 450 core

in vec3 normal;
out vec4 color;

uniform vec4 _color;

void main(void) {
	color = _color;
}
