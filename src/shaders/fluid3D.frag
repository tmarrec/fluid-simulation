#version 460 core

out vec4 FragColor;

in VS_OUT
{
    vec3 FragPos;
    vec3 Normal;
    vec3 TexCoord;
} fs_in;

uniform sampler3D densityTex;
uniform vec3 eyePos;

void main()
{
    FragColor = vec4(1, 1, 1, 1.0f);
}
