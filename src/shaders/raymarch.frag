#version 460 core

out vec4 FragColor;

uniform sampler3D u_densityTex;
uniform vec2 u_resolution;  // Width and height of the shader
uniform vec3 u_eyePos; // Camera positi1n
uniform vec3 u_eyeFront; // Camera front
uniform float u_eyeFOV; // Camera field of view
uniform float u_absorption;
uniform vec3 u_lightIntensity;

#define MAX_STEPS 64
#define MAX_DIST 64.0f
#define SURFACE_DIST 0.01f

const int numSamples = 128;
const float cubeSize = 3.0f;
const float maxDist = sqrt(3.0f);
const float scale = maxDist / float(numSamples);
const int numLightSamples = 64;
const float lscale = maxDist / float(numLightSamples);
const vec3 lightPos = vec3(2,4,1);

float getCubeDist(vec3 p)
{
	vec3 s = vec3(cubeSize);
    vec3 bp = p;
    p = abs(p)-s;
	float cube = length(max(p, 0.0))+min(max(p.x, max(p.y, p.z)), 0.);
	return cube;
}
 
vec4 rayMarching(vec3 ro, vec3 rd) 
{
	float dO = 0.0f; // Distance origin
	for(int i = 0; i < MAX_STEPS; ++i)
	{
		vec3 p = ro + rd * dO;
		float ds = getCubeDist(p); // ds is Distance Scene
		dO += ds;
        if (dO > MAX_DIST)
            discard;
		if (ds < SURFACE_DIST) // If on cube surface
        {
            vec3 pos = (p+vec3(cubeSize))/(cubeSize*2); // Tex pos
            vec3 eyeDir = normalize(pos-u_eyePos)*scale;

            float T = 0.8f; // Transmittance
            vec3 Lo = vec3(0.0f); // In-scattered radiance

            for (int s = 0; s < numSamples; ++s)
            {
                float density = texture(u_densityTex, pos).x;

                if (density > 0)
                {
                    T *= 1.0f-density*scale*u_absorption;
                    if (T <= 0.01f)
                        break;

                    vec3 lightDir = normalize(lightPos-pos)*lscale;

                    float Tl = 1.0f;
                    vec3 lpos = pos + lightDir;

                    for (int sl = 0; sl < numLightSamples; ++sl)
                    {
                        float ld = texture(u_densityTex, lpos).x;
                        Tl *= 1.0f-u_absorption*ld*lscale;

                        if (Tl <= 0.01f)
                            break;

                        lpos += lightDir;
                    }

                    vec3 Li = u_lightIntensity*Tl*u_absorption;
                    Lo += Li*T*density*scale;
                }

                pos += eyeDir; // Going deeper inside the cube
            }
            return vec4(Lo, 1.0f-T);
        }
    }
    discard;
}

vec3 R(vec2 uv, vec3 p, vec3 l, float z)
{
    vec3 f = normalize(l-p);
    vec3 r = normalize(cross(vec3(0,1,0), f));
    vec3 u = cross(f, r);
    vec3 c = p+f*z;
    vec3 i = c + -uv.x*r + uv.y*u;
    vec3 d = normalize(i-p);
    return d;
}
 
void main()
{
	vec2 uv = (gl_FragCoord.xy-0.5f*u_resolution.xy)/u_resolution.y;
	vec3 ro = u_eyePos; // Ray origin
    //TODO fix FOV under here
    vec3 rd = R(uv, ro, u_eyePos+u_eyeFront, radians(u_eyeFOV-10)); // Ray direction
    FragColor = rayMarching(ro, rd);
}

