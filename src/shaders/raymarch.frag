#version 460 core

in vec2 TexCoords;

out vec4 FragColor;

uniform sampler2D screenTexture;

/*
uniform sampler3D densityTex;
uniform vec3 lightPos;
uniform vec3 lightIntensity;
uniform vec3 eyePos;
uniform float absorption;
*/

uniform vec2 u_resolution;  // Width and height of the shader
uniform float u_time;  // Time elapsed
uniform vec3 u_eyePos; // Camera position

#define MAX_STEPS 100
#define MAX_DIST 100.0f
#define SURFACE_DIST 0.01f

mat2 Rot(float a)
{
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

float GetDist(vec3 p)
{
	vec4 s = vec4(0,2,8,2); //Sphere. xyz is position w is radius
	float sphereDist = length(p-s.xyz) - s.w;
	float planeDist = p.y;
	float d = min(sphereDist, sphereDist);
 
	return d;
}
 
float RayMarch(vec3 ro, vec3 rd) 
{
	float dO = 0.0f; //Distane Origin
	for(int i = 0; i < MAX_STEPS; ++i)
	{
		vec3 p = ro + rd * dO;
		float ds = GetDist(p); // ds is Distance Scene
		dO += ds;
        if (dO > MAX_DIST)
            discard;
		if (ds < SURFACE_DIST)
            break;
    }
	return dO;
}

vec3 GetNormal(vec3 p)
{ 
    float d = GetDist(p); // Distance
    vec2 e = vec2(.01,0); // Epsilon
    vec3 n = d - vec3(
    GetDist(p-e.xyy),  
    GetDist(p-e.yxy),
    GetDist(p-e.yyx));
   
    return normalize(n);
}

float GetLight(vec3 p)
{ 
    // Light (directional diffuse)
    vec3 lightPos = vec3(5.*sin(u_time),5.,5.0*cos(u_time)); // Light Position
    vec3 l = normalize(lightPos-p); // Light Vector
    vec3 n = GetNormal(p); // Normal Vector
   
    float dif = dot(n,l); // Diffuse light
    dif = clamp(dif,0.,1.); // Clamp so it doesnt go below 0
 
    return dif;
}
 
void main()
{
	vec2 uv = (gl_FragCoord.xy-0.5f*u_resolution.xy)/u_resolution.y;
	vec3 ro = u_eyePos; // Ray Origin/ Camera
	vec3 rd = normalize(vec3(uv.x, uv.y, 1.0f));

	float d = RayMarch(ro, rd); // Distance

    vec3 p = ro + rd * d;
    float dif = GetLight(p);
    d *= 0.2f;

    vec3 color = vec3(dif);

	FragColor = vec4(color, 0.9f);
}

/*
void main()
{
	// Diagonal of the cube
	const float maxDist = sqrt(3.0f);

	const int numSamples = 32;
	const float scale = maxDist/float(numSamples);

	const int numLightSamples = 8;
	const float lscale = maxDist/float(numLightSamples);

	vec3 pos = vec3(TexCoords.xy, 1);
	vec3 eyeDir = normalize(pos-eyePos)*scale;

	float T = 1.0f;
	vec3 Lo = vec3(0.0f);

	for (int i = 0; i < numSamples; ++i)
	{
		float density = texture(densityTex, pos).x;
		if (density > 0.0f)
		{
			T *= 1.0-density*scale*absorption;
			if (T <= 0.01)
			{
				break;
			}

			vec3 lightDir = normalize(lightPos-pos)*lscale;

			float Tl = 1.0f;
			vec3 lpos = pos + lightDir;

			for (int s = 0; s < numLightSamples; ++s)
			{
				float ld = texture(densityTex, lpos).x;
				Tl *= 1.0f-absorption*lscale*ld;

				if (Tl <= 0.01f)
				{
					break;
				}

				lpos += lightDir;
			}

			vec3 Li = lightIntensity*Tl;
			Lo += Li*T*density*scale;
		}

		pos += eyeDir;
	}

	FragColor.xyz = Lo;
	FragColor.w = 1.0f-T;
}
*/

