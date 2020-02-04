#version 450 core

in vec3 normal;
out vec4 color;

uniform vec4 time;

void main(void) {
	//color = vec4(1.0f, 0.2f, 0.2f, 1.0f);
	color = time;
}
