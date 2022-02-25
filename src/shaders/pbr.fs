#version 330 core
out vec4 FragColor;

in vec3 FragPos;
in vec2 TexCoords;
in vec3 Normal;
in mat3 TBN;

const float PI = 3.14159265359;

struct Material {
    vec3 material_ambient;
    vec3 material_diffuse;
    vec3 material_specular;
    vec3 material_emmissive;
    float material_shininess;
    float material_refracti;
}

uniform vec3 material_ambient;
uniform vec3 material_diffuse;
uniform vec3 material_specular;
uniform float material_shininess;
uniform float material_refracti;

uniform bool use_texture_ambient;
uniform bool use_texture_diffuse;
uniform bool use_texture_specular;
uniform bool use_texture_normal;
uniform bool use_texture_height;

uniform sampler2D texture_ambient;
uniform sampler2D texture_diffuse;
uniform sampler2D texture_specular;
uniform sampler2D texture_normal;
uniform sampler2D texture_height;

uniform vec3 lightPos;
uniform vec3 viewPos;

vec2 ParallaxMapping(vec2 texCoords, vec3 viewDir)
{
    float height_scale = 0.1;
    // number of depth layers
    const float minLayers = 8;
    const float maxLayers = 32;
    float numLayers = mix(maxLayers, minLayers, abs(dot(vec3(0.0, 0.0, 1.0), viewDir)));  
    // calculate the size of each layer
    float layerDepth = 1.0 / numLayers;
    // depth of current layer
    float currentLayerDepth = 0.0;
    // the amount to shift the texture coordinates per layer (from vector P)
    vec2 P = viewDir.xy / viewDir.z * height_scale;
    vec2 deltaTexCoords = P / numLayers;
  
    // get initial values
    vec2 currentTexCoords = texCoords;
    float currentDepthMapValue = texture(texture_height, currentTexCoords).r;
      
    while (currentLayerDepth < currentDepthMapValue) {
        // shift texture coordinates along direction of P
        currentTexCoords -= deltaTexCoords;
        // get depthmap value at current texture coordinates
        currentDepthMapValue = texture(texture_height, currentTexCoords).r;  
        // get depth of next layer
        currentLayerDepth += layerDepth;  
    }
    
    // -- parallax occlusion mapping interpolation from here on
    // get texture coordinates before collision (reverse operations)
    vec2 prevTexCoords = currentTexCoords + deltaTexCoords;

    // get depth after and before collision for linear interpolation
    float afterDepth  = currentDepthMapValue - currentLayerDepth;
    float beforeDepth = texture(texture_height, prevTexCoords).r - currentLayerDepth + layerDepth;
 
    // interpolation of texture coordinates
    float weight = afterDepth / (afterDepth - beforeDepth);
    vec2 finalTexCoords = prevTexCoords * weight + currentTexCoords * (1.0 - weight);

    return finalTexCoords;
}

// -----------------------------------------------
float DistributionGGX(vec3 N, vec3 H, float roughness)
{
    float a = roughness * roughness;
    float a2 = a * a;
    float NdotH = max(dot(N, H), 0.0);
    float NdotH2 = NdotH*NdotH;

    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;

    return a2 / denom;
}
float GeometrySmithGGX(vec3 N, vec3 V, vec3 L, float roughness)
{
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;

    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx2 = NdotV / (NdotV * (1.0 - k) + k);
    float ggx1 = NdotL / (NdotL * (1.0 - k) + k);

    return ggx1 * ggx2;
}
vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}
vec4 BRDF(vec3 viewDir, vec3 lightDir, vec3 normal, vec3 diffuse_color, float roughness, float distance)
{
    vec3 F0 = vec3(pow((material_refracti - 1) / (material_refracti + 1), 2));
    float attenuation = 1.0 / (distance * distance);
    vec3 halfwayDir = normalize(lightDir + viewDir);

    // Cook-Torrance BRDF
    float D = DistributionGGX(normal, halfwayDir, roughness);   
    float G = GeometrySmithGGX(normal, viewDir, lightDir, roughness);      
    vec3 F  = fresnelSchlick(clamp(dot(halfwayDir, viewDir), 0.0, 1.0), F0);
           
    vec3 numerator = D * G * F; 
    float denominator = 4.0 * max(dot(normal, viewDir), 0.0) * max(dot(normal, lightDir), 0.0) + 0.0001; // + 0.0001 to prevent divide by zero
    vec3 specular = numerator / denominator;
        
    // kS is equal to Fresnel
    vec3 kS = F;
    vec3 kD = vec3(1.0) - kS; 

    // scale light by NdotL
    float NdotL = max(dot(normal, lightDir), 0.0);
    vec3 Lo = (kD * diffuse_color / PI + specular) * attenuation * NdotL;
    vec3 ambient = vec3(0.03) * diffuse_color;

    vec3 color = ambient + Lo;
    color = color / (color + vec3(1.0));
    color = pow(color, vec3(1.0/2.2));

    return vec4(color, 1.0);
}

void main()
{
    vec3 view_pos = viewPos;
    vec3 light_pos = lightPos;
    vec3 frag_pos = FragPos;
    vec3 ambient_color = material_ambient;
    vec3 diffuse_color = material_diffuse;
    vec3 specular_color = material_specular;
    vec2 texCoords = TexCoords;
    vec3 normal = Normal;

    // check if need to use tangent space
    if (use_texture_height || use_texture_normal) {
        normal = vec3(0, 0, 1);
        light_pos = TBN * lightPos;
        view_pos  = TBN * viewPos;
        frag_pos  = TBN * FragPos;
    }
    vec3 viewDir = normalize(view_pos - frag_pos);
    vec3 lightDir = normalize(light_pos - frag_pos);
    float distance = length(light_pos - frag_pos);

    // check if use textures
    if (use_texture_height) {
        texCoords = ParallaxMapping(TexCoords, viewDir);
        if (texCoords.x > 1.0 || texCoords.y > 1.0 || texCoords.x < 0.0 || texCoords.y < 0.0)
            discard;
    }
    if (use_texture_normal) {
        normal = texture(texture_normal, texCoords).rgb;
        normal = normalize(normal * 2.0 - 1.0);
    }
    if (use_texture_ambient) {
        ambient_color = ambient_color * texture(texture_ambient, texCoords).rgb;
    }
    if (use_texture_diffuse) {
        diffuse_color = diffuse_color * texture(texture_diffuse, texCoords).rgb;
    }
    if (use_texture_specular) {
        specular_color = specular_color * texture(texture_specular, texCoords).rgb;
    }

    FragColor = BRDF(viewDir, lightDir, normal, diffuse_color, (1.0 - specular_color.r) * (1.0 - material_shininess / 1000), distance);
    // FragColor = vec4(vec3((1.0 - specular_color.r) * (1.0 - material_shininess / 1000)), 1.0);
} 