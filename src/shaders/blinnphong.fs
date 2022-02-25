#version 330 core
out vec4 FragColor;

in vec3 FragPos;
in vec2 TexCoords;
in vec3 Normal;
in mat3 TBN;

struct Material {
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
    vec3 emmissive;
    float shininess;
    float metallic;
    float refracti;
}

uniform Material material;

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

vec4 BlinnPhong(vec3 viewDir, vec3 lightDir, vec3 normal, vec3 ambient_color, vec3 diffuse_color, vec3 specular_color, float shininess, float distance)
{
    // ambient
    vec3 ambient = 0.1 * ambient_color * diffuse_color;
    // diffuse
    float diff = max(dot(normal, lightDir), 0.0);
    vec3 diffuse = diffuse_color * diff;  
    // specular
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(normal, halfwayDir), 0.0), shininess);
    vec3 specular = specular_color * spec;
    
    // attenuation
    float attenuation = 1.0 / (1 + 0.09 * distance + 0.032 * (distance * distance));
    ambient *= attenuation;
    diffuse *= attenuation;
    specular *= attenuation;

    return vec4(ambient + diffuse + specular, 1.0);
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

    // ambient
    vec3 ambient = 0.1 * ambient_color * diffuse_color;
    // diffuse
    float diff = max(dot(normal, lightDir), 0.0);
    vec3 diffuse = diffuse_color * diff;  
    // specular
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(normal, halfwayDir), 0.0), material_shininess);
    vec3 specular = specular_color * spec;
    
    // attenuation
    float attenuation = 1.0 / (1 + 0.09 * distance + 0.032 * (distance * distance));
    ambient *= attenuation;
    diffuse *= attenuation;
    specular *= attenuation;

    FragColor = vec4(ambient + diffuse + specular, 1.0);
} 