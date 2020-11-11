#version 460 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out VS_OUT
{
	vec3 Normal;
} vs_out;

void main(void)
{
	vs_out.Normal = aNormal;
	gl_Position = vec4(aPos, 1.0);
}
