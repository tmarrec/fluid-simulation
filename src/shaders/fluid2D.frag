#version 460 core

out vec4 FragColor;

in VS_OUT
{
    vec3 FragPos;
    vec3 Normal;
    vec2 TexCoord;
} fs_in;

uniform sampler2D tex0;

void main()
{
    //FragColor = vec3(fs_in.TexCoord, 0, 0.0f);
    FragColor = texture(tex0, fs_in.TexCoord);
    //FragColor = vec4(0, 0, 0, 1.0f);
}
