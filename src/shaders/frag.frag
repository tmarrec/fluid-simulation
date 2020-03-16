#version 410 core

#define N_MAX_LIGHT 128

struct Point_Light {
    vec3 position;
    vec3 color;

	float intensity;
	float constant;
	float linear;
	float quadratic;
};


in vec3 normal;
in vec3 pos;
out vec4 color;

uniform vec3 _object_color;
uniform vec3 _view_pos;
uniform int _light_nb;

uniform Point_Light _point_lights[N_MAX_LIGHT];

// Functions
vec3 calc_point_light(Point_Light light, vec3 normal, vec3 frag_pos, vec3 view_dir);

void main(void) {
	vec3 norm = normalize(normal);
	vec3 view_dir = normalize(_view_pos - pos);

	vec3 result = vec3(0,0,0);

	// Limits lights number
	int light_nb = _light_nb;
	if (_light_nb > N_MAX_LIGHT) {
		light_nb = N_MAX_LIGHT;
	}

	// Compute light impact
	for (int i = 0; i < light_nb; i++) {
		result += calc_point_light(_point_lights[i], norm, pos, view_dir);
	}

	color = vec4(result * _object_color, 1.0);
}

vec3 calc_point_light(Point_Light light, vec3 normal, vec3 frag_pos, vec3 view_dir) {
	// Attenuation
	float distance    = length(light.position - frag_pos);
    float attenuation = 1.0 / (light.constant + light.linear * distance + light.quadratic * (distance * distance));    

	// Ambient
	float ambient_amp = 0.1;
	vec3 ambient = ambient_amp * light.color * light.intensity;

	// Diffuse
	vec3 light_dir = normalize(light.position - frag_pos);
	float diff = max(dot(normal, light_dir), 0.0);
	vec3 diffuse = diff * light.color * light.intensity;

	// Specular
	float specular_strength = 0.5;
	vec3 reflect_dir = reflect(-light_dir, normal);
	float spec = pow(max(dot(view_dir, reflect_dir), 0.0), 8);
	vec3 specular = specular_strength * spec * light.color * light.intensity;

	ambient *= attenuation;
	diffuse *= attenuation;
	specular *= attenuation;
	

	return (ambient + diffuse + specular);
}
