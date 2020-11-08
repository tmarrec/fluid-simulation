#version 460 core

#define N_MAX_LIGHT 16

struct Point_Light
{
	vec3 position;
	vec3 color;

	float intensity;
	float constant;
	float linear;
	float quadratic;
};

in VS_OUT
{
	vec3 FragPos;
	vec3 Normal;
} fs_in;

out vec4 FragColor;

uniform vec3 _object_color;
uniform vec3 _view_pos;
uniform int _light_nb;

uniform float farPlane;
uniform samplerCube shadowMaps[N_MAX_LIGHT];

uniform Point_Light _point_lights[N_MAX_LIGHT];


float ShadowMap(vec3 normal, vec3 lightDir, vec3 lightPos, samplerCube shadowMap)
{
	// Perspective divide

	vec3 fragToLight = fs_in.FragPos - lightPos;
	float closestDepth = texture(shadowMap, fragToLight).r;
	closestDepth *= farPlane;
	float currentDepth = length(fragToLight);
	float acneeBias = max(0.05*(1.0-dot(normal, lightDir)), 0.005);
	float shadow = (currentDepth - acneeBias) > closestDepth ? 1.0 : 0.0;
	return shadow;
}

vec3 calc_point_light(int lightInd, Point_Light light, vec3 normal, vec3 frag_pos, vec3 view_dir)
{
	// Attenuation
	float distance	= length(light.position - frag_pos);
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

	// Shadow Map
	float shadow = ShadowMap(normal, light_dir, light.position, shadowMaps[lightInd]);

	ambient *= attenuation;
	diffuse *= attenuation;
	specular *= attenuation;

	return (ambient + (1.0 - shadow) * (diffuse + specular));
}

void main(void)
{
	vec3 norm = normalize(fs_in.Normal);
	vec3 view_dir = normalize(_view_pos - fs_in.FragPos);

	vec3 result = vec3(0,0,0);

	// Limits lights number
	int light_nb = _light_nb;
	if (_light_nb > N_MAX_LIGHT)
	{
		light_nb = N_MAX_LIGHT;
	}

	// Compute light impact
	for (int i = 0; i < light_nb; i++)
	{
		result += calc_point_light(i, _point_lights[i], norm, fs_in.FragPos, view_dir);
	}

	FragColor = vec4(result * _object_color, 1.0);
}
