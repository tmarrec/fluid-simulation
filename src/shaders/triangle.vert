#version 450 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 inormal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec3 normal;
out vec3 pos;

void main(void) {
	normal = inormal;
	pos = position;
	gl_Position = projection * view * model * vec4(position, 1.0f);
}
