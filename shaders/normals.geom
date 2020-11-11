#version 460

layout (triangles) in;
layout (line_strip, max_vertices=6) out;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

in VS_OUT
{
	vec3 Normal;
} vs_out[];

void main()
{
	int i;
    for(i = 0; i < gl_in.length(); i++)
    {
		vec3 P = gl_in[i].gl_Position.xyz;
		vec3 N = vs_out[i].Normal;

        gl_Position = projection * view * model * vec4(P, 1.0);
        EmitVertex();

		float normalLength = 0.5;
        gl_Position = projection * view * model * vec4(P+N*normalLength, 1.0);
        EmitVertex();

    	EndPrimitive();
    }
}
