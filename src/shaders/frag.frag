#version 460 core

#define N_MAX_LIGHT 24

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
	vec4 FragPosLightSpace[N_MAX_LIGHT];
} fs_in;

out vec4 FragColor;

uniform vec3 _object_color;
uniform vec3 _view_pos;
uniform int _light_nb;

uniform int shadowMapNb;
uniform sampler2D shadowMaps[N_MAX_LIGHT];

uniform Point_Light _point_lights[N_MAX_LIGHT];

// Functions
vec3 calc_point_light(int lightInd, Point_Light light, vec3 normal, vec3 frag_pos, vec3 view_dir, vec4 fragPosLightSpace);

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
		result += calc_point_light(i, _point_lights[i], norm, fs_in.FragPos, view_dir, fs_in.FragPosLightSpace[i]);
	}

	FragColor = vec4(result * _object_color, 1.0);
}

float ShadowMap(vec4 fragPosLightSpace, vec3 normal, vec3 lightDir, sampler2D shadowMap)
{
	// Perspective divide
	vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;
	projCoords = projCoords * 0.5 + 0.5; // [-1;1] -> [0;1]
	if (projCoords.z > 1.0)
	{
		return 0.0f;
	}

	float closestDepth = texture(shadowMap, projCoords.xy).r;
	float currentDepth = projCoords.z;
	float acneeBias = max(0.000005*(1.0-dot(normal, lightDir)), 0.0000005);

	float shadow = 0.0;
	vec2 texelSize = 1.0 / textureSize(shadowMap, 0);
	for (int x = -1; x <= 1; ++x)
	{
		for (int y = -1; y <= 1; ++y)
		{
			float pcfDepth = texture(shadowMap, projCoords.xy + vec2(x,y) * texelSize).r;
			shadow += (currentDepth - acneeBias) > pcfDepth ? 1.0 : 0.0;
		}
	}

	return shadow /= 9.0;
}

vec3 calc_point_light(int lightInd, Point_Light light, vec3 normal, vec3 frag_pos, vec3 view_dir, vec4 fragPosLightSpace)
{
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

	// Shadow Map
	float shadow = ShadowMap(fragPosLightSpace, normal, light_dir, shadowMaps[lightInd]);

	ambient *= attenuation;
	diffuse *= attenuation;
	specular *= attenuation;

	return (ambient + (1.0 - shadow) * (diffuse + specular));
}
