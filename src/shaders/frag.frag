#version 450 core

in vec3 normal;
in vec3 pos;
out vec4 color;

uniform vec3 _object_color;
uniform vec3 _view_pos;
uniform vec3 _light_color;
uniform vec3 _light_pos;

void main(void) {
	// ambient
	float ambient_amp = 0.1;
	vec3 ambient = ambient_amp * _light_color;

	// diffuse
	vec3 norm = normalize(normal);
	vec3 lightDir = normalize(_light_pos - pos);
	float diff = max(dot(norm, lightDir), 0.0);
	vec3 diffuse = diff * _light_color;

	// specular
	float specular_strength = 0.5;
	vec3 view_dir = normalize(_view_pos - pos);
	vec3 reflect_dir = reflect(-lightDir, norm);
	float spec = pow(max(dot(view_dir, reflect_dir), 0.0), 16);
	vec3 specular = specular_strength * spec * _light_color;

	vec3 result = (ambient + diffuse + specular) * _object_color;
	color = vec4(result, 1.0);
}
