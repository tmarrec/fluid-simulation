#version 450 core

in vec3 normal;
in vec3 pos;
out vec4 color;

uniform vec3 _color;

void main(void) {
	// ambient
	float ambient_amp = 0.1;
	vec3 ambient = ambient_amp * vec3(1,1,1); // <- light color

	// diffuse
	vec3 norm = normalize(normal);
	vec3 lightDir = normalize(vec3(0,0,10) - pos); // <- light pos
	float diff = max(dot(norm, lightDir), 0.0);
	vec3 diffuse = diff * vec3(1,1,1); // <- light color

	vec3 result = (ambient + diffuse) * _color;
	color = vec4(result, 1.0);
}
