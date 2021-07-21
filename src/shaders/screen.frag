#version 460 core

in vec2 TexCoords;

out vec4 FragColor;

uniform sampler2D screenTexture;

void main()
{
	vec3 col = texture(screenTexture, TexCoords).rgb;
	FragColor = vec4(col, 1.0f);
}
