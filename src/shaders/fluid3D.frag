#version 460 core

out vec4 FragColor;

in VS_OUT
{
    vec3 FragPos;
    vec3 Normal;
    vec3 TexCoord;
} fs_in;

uniform sampler3D densityTex;
uniform vec3 eyePos;

void main()
{
    //FragColor = texture(densityTex, fs_in.TexCoord);
    FragColor = vec4(1, 1, 1, 1.0f);
    //
    //

    /*
    // Diagonal of the cube
    const float maxDist = sqrt(3.0f);

    const int numSamples = 128;
    const float scale = maxDist/float(numSamples);

    const int numLightSamples = 32;
    const float lscale = maxDist/float(numLightSamples);

    vec3 pos = fs_in.TexCoord.xyz;
    vec3 eyeDir = normalize(pos-eyePos)*scale;

    float T = 1.0f;
    vec3 Lo = vec3(0.0f);


    // TODO TEMP
    float absorption = 1.0f;
    vec3 lightPos = vec3(0,0,0);
    vec3 lightIntensity = vec3(1,1,1);

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
    */
}
