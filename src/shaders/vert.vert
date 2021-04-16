#version 460 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 3) in vec3 aTexCoord;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out VS_OUT
{
    vec3 FragPos;
    vec3 Normal;
    vec3 TexCoord;
} vs_out;

void main()
{
    vs_out.FragPos = vec3(model*vec4(aPos, 1.0));
    vs_out.Normal = mat3(transpose(inverse(model))) * aNormal;
    vs_out.TexCoord = aTexCoord;
    gl_Position = projection * view * vec4(vs_out.FragPos, 1.0);
}
