#version 460 core

#define N_MAX_LIGHT 24

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform mat4 lightSpaceMatrix[N_MAX_LIGHT];
uniform int shadowMapNb;

out VS_OUT
{
	vec3 FragPos;
	vec3 Normal;
	vec4 FragPosLightSpace[N_MAX_LIGHT];
} vs_out;

void main(void) {
	vs_out.FragPos = vec3(model * vec4(aPos, 1.0));
	vs_out.Normal = mat3(transpose(inverse(model))) * aNormal;
	for (int i = 0; i < shadowMapNb; i++)
	{
		vs_out.FragPosLightSpace[i] = lightSpaceMatrix[i] * vec4(vs_out.FragPos, 1.0);
	}
	gl_Position = projection * view * vec4(vs_out.FragPos, 1.0);
}
