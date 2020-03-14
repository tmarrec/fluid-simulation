#version 450 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 i_normal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec3 normal;
out vec3 pos;
out vec3 pos_error;

void main(void) {
	normal = mat3(transpose(inverse(model)))*i_normal;
	pos = vec3(model*vec4(position, 1.0));
	pos_error = position;
	gl_Position = projection * view * model * vec4(position, 1.0f);
}
