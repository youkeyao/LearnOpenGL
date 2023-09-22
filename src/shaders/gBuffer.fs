#version 330 core
layout (location = 0) out vec3 fragPos;
layout (location = 1) out vec3 fragNorm;
layout (location = 2) out vec3 fragAlbedo;
layout (location = 3) out vec3 motionVecRoughness;

in vec3 v_Position;
in vec3 v_Normal;
in vec2 v_Texcoord;
in float v_MatID;
in vec4 nowVertex;
in vec4 prevVertex;

uniform sampler2DArray texturesArray;
uniform samplerBuffer materials;

struct Material {
    vec3 ambient;       // texture: (-id-2, width, height) color: (r, g, b)
    vec3 diffuse;
    vec3 specular;
    vec3 emissive;
    float shininess;
    float roughness;
    float metallic;
    float refracti;
    float opacity;
    vec3 transmission;
    float anisotropy;
};

vec3 getColor(int pos, vec2 texCoords)
{
    vec3 t = texelFetch(materials, pos).xyz;
    if (t.x < 0) {
        return texture2DArray(texturesArray, vec3(texCoords * t.yz, - t.x - 2)).rgb;
    }
    else {
        return t;
    }
}

void main()
{
    fragPos = v_Position;
    fragNorm = v_Normal;
    fragAlbedo = getColor(int(v_MatID) * 10 + 1, v_Texcoord);
    vec2 motionVec = (nowVertex.rg/nowVertex.w - prevVertex.rg/prevVertex.w) / 2;
    motionVecRoughness = vec3(motionVec, 1 - texelFetch(materials, int(v_MatID) * 10 + 4).x / 1000);
}