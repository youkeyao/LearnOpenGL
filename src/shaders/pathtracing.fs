#version 330 core

#define SAMPLENUM 1

out vec4 FragColor;

uniform float screenWidth;
uniform float screenHeight;
uniform mat4 projection;
uniform mat4 view;

uniform uint counter;

uniform sampler2DArray texturesArray;
uniform samplerBuffer triangles;
uniform int nTriangles;

uniform sampler2D lastFrame;

uniform vec3 viewPos;

const float PI = 3.14159265359;

struct Material {
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
    vec3 emissive;
    vec3 shininess;
    vec3 metallic;
    vec3 refracti;
    vec3 opacity;
    vec3 transmission;
    vec3 anisotropy;
};

struct Triangle {
    vec3 p[3];
    vec3 n[3];
    vec2 texCoords[3];
    mat3 TBN[3];
    vec3 normal;
    vec3 height;
};

struct Ray {
    vec3 startPoint;
    vec3 direction;
};

struct HitResult {
    bool isHit;
    bool isInside;
    float distance;
    vec3 hitPoint;
    vec3 normal;
    vec2 texCoords;
    vec3 viewDir;
    mat3 TBN;
    Material material;
};

Triangle getTriangle(int i)
{
    int offset = i * 29;
    Triangle t;

    t.p[0] = texelFetch(triangles, offset + 0).xyz;
    t.p[1] = texelFetch(triangles, offset + 1).xyz;
    t.p[2] = texelFetch(triangles, offset + 2).xyz;

    t.n[0] = texelFetch(triangles, offset + 3).xyz;
    t.n[1] = texelFetch(triangles, offset + 4).xyz;
    t.n[2] = texelFetch(triangles, offset + 5).xyz;

    t.texCoords[0] = texelFetch(triangles, offset + 6).xy;
    t.texCoords[1] = vec2(texelFetch(triangles, offset + 6).z, texelFetch(triangles, offset + 7).x);
    t.texCoords[2] = texelFetch(triangles, offset + 7).yz;

    t.TBN[0] = mat3(texelFetch(triangles, offset + 8).xyz, texelFetch(triangles, offset + 9).xyz, texelFetch(triangles, offset + 10).xyz);
    t.TBN[1] = mat3(texelFetch(triangles, offset + 11).xyz, texelFetch(triangles, offset + 12).xyz, texelFetch(triangles, offset + 13).xyz);
    t.TBN[2] = mat3(texelFetch(triangles, offset + 14).xyz, texelFetch(triangles, offset + 15).xyz, texelFetch(triangles, offset + 16).xyz);

    t.normal = texelFetch(triangles, offset + 17).xyz;
    t.height = texelFetch(triangles, offset + 18).xyz;

    return t;
}

Material getMaterial(int i, vec2 texCoords)
{
    int offset = i * 29;
    Material m;
    vec3 tmp;

    tmp = texelFetch(triangles, offset + 19).xyz;
    m.ambient = texture2DArray(texturesArray, vec3(texCoords * tmp.yz, tmp.x)).rgb;
    tmp = texelFetch(triangles, offset + 20).xyz;
    m.diffuse = texture2DArray(texturesArray, vec3(texCoords * tmp.yz, tmp.x)).rgb;
    tmp = texelFetch(triangles, offset + 21).xyz;
    m.specular = texture2DArray(texturesArray, vec3(texCoords * tmp.yz, tmp.x)).rgb;
    tmp = texelFetch(triangles, offset + 22).xyz;
    m.emissive = texture2DArray(texturesArray, vec3(texCoords * tmp.yz, tmp.x)).rgb;
    tmp = texelFetch(triangles, offset + 23).xyz;
    m.shininess = texture2DArray(texturesArray, vec3(texCoords * tmp.yz, tmp.x)).rgb;
    tmp = texelFetch(triangles, offset + 24).xyz;
    m.metallic = texture2DArray(texturesArray, vec3(texCoords * tmp.yz, tmp.x)).rgb;
    tmp = texelFetch(triangles, offset + 25).xyz;
    m.refracti = texture2DArray(texturesArray, vec3(texCoords * tmp.yz, tmp.x)).rgb;
    tmp = texelFetch(triangles, offset + 26).xyz;
    m.opacity = texture2DArray(texturesArray, vec3(texCoords * tmp.yz, tmp.x)).rgb;
    tmp = texelFetch(triangles, offset + 27).xyz;
    m.transmission = texture2DArray(texturesArray, vec3(texCoords * tmp.yz, tmp.x)).rgb;
    tmp = texelFetch(triangles, offset + 28).xyz;
    m.anisotropy = texture2DArray(texturesArray, vec3(texCoords * tmp.yz, tmp.x)).rgb;

    return m;
}

vec3 calcLerp(vec3 p1, vec3 p2, vec3 p3, vec3 P)
{
    // float a = (-(P.x-p2.x)*(p3.y-p2.y) + (P.y-p2.y)*(p3.x-p2.x)) / (-(p1.x-p2.x-0.00005)*(p3.y-p2.y+0.00005) + (p1.y-p2.y+0.00005)*(p3.x-p2.x+0.00005));
    // float b  = (-(P.x-p3.x)*(p1.y-p3.y) + (P.y-p3.y)*(p1.x-p3.x)) / (-(p2.x-p3.x-0.00005)*(p1.y-p3.y+0.00005) + (p2.y-p3.y+0.00005)*(p1.x-p3.x+0.00005));
    // float c  = 1.0 - a - b;

    // return vec3(a, b, c);
    float a1, a2, b1, b2, c1, c2;

    if (abs((p1.x - p3.x) * (p2.y - p3.y) - (p1.y - p3.y) * (p2.x - p3.x)) > 0.00001f) {
        a1 = p1.x - p3.x;
        a2 = p1.y - p3.y;
        b1 = p2.x - p3.x;
        b2 = p2.y - p3.y;
        c1 = P.x - p3.x;
        c2 = P.y - p3.y;
    }
    else if (abs((p1.x - p3.x) * (p2.z - p3.z) - (p1.z - p3.z) * (p2.x - p3.x)) > 0.00001f) {
        a1 = p1.x - p3.x;
        a2 = p1.z - p3.z;
        b1 = p2.x - p3.x;
        b2 = p2.z - p3.z;
        c1 = P.x - p3.x;
        c2 = P.z - p3.z;
    }
    else {
        a1 = p1.y - p3.y;
        a2 = p1.z - p3.z;
        b1 = p2.y - p3.y;
        b2 = p2.z - p3.z;
        c1 = P.y - p3.y;
        c2 = P.z - p3.z;
    }

    float d = a1 * b2 - a2 * b1;
    float a = (c1 * b2 - c2 * b1) / d;
    float b = (a1 * c2 - a2 * c1) / d;
    float c = 1.0 - a - b;

    return vec3(a, b, c);
}

HitResult hitTriangle(Triangle triangle, Ray ray)
{
    HitResult res;
    res.distance = 0;
    res.isHit = false;
    res.isInside = false;

    vec3 p1 = triangle.p[0];
    vec3 p2 = triangle.p[1];
    vec3 p3 = triangle.p[2];

    vec3 S = ray.startPoint;
    vec3 d = ray.direction;
    vec3 N = normalize(cross(p2-p1, p3-p1));

    if (dot(N, d) > 0.0f) {
        N = -N;   
        res.isInside = true;
    }

    if (abs(dot(N, d)) < 0.00001f) return res;

    float t = (dot(N, p1) - dot(S, N)) / dot(d, N);
    if (t < 0.0005f) return res;

    vec3 P = S + d * t;

    vec3 c1 = cross(p2 - p1, P - p1);
    vec3 c2 = cross(p3 - p2, P - p2);
    vec3 c3 = cross(p1 - p3, P - p3);
    bool r1 = (dot(c1, N) > -0.01f && dot(c2, N) > -0.01f && dot(c3, N) > -0.01f);
    bool r2 = (dot(c1, N) < 0.01f && dot(c2, N) < 0.01f && dot(c3, N) < 0.01f);

    if (r1 || r2) {
        res.isHit = true;
        res.hitPoint = P;
        res.distance = t;
        res.normal = N;
        res.viewDir = d;

        vec3 lerp = calcLerp(p1, p2, p3, P);
        res.texCoords = lerp.x * triangle.texCoords[0] + lerp.y * triangle.texCoords[1] + lerp.z * triangle.texCoords[2];
        res.TBN = lerp.x * triangle.TBN[0] + lerp.y * triangle.TBN[1] + lerp.z * triangle.TBN[2];

        if (triangle.normal.x >= 0) {
            res.normal = texture2DArray(texturesArray, vec3(res.texCoords * triangle.normal.yz, triangle.normal.x)).rgb;
            res.normal = normalize(res.normal * 2.0 - 1.0);
            res.normal = res.TBN * res.normal;
        }
        else {
            res.normal = lerp.x * triangle.n[0] + lerp.y * triangle.n[1] + lerp.z * triangle.n[2];
        }
        res.normal = normalize(res.normal);
        res.normal = (res.isInside) ? (-res.normal) : (res.normal);
    }

    return res;
}

HitResult hitArray(Ray ray, int l, int r)
{
    HitResult res;
    res.isHit = false;
    for (int i = l; i <= r; i ++) {
        Triangle triangle = getTriangle(i);
        HitResult r = hitTriangle(triangle, ray);
        if (r.isHit && (!res.isHit || r.distance < res.distance)) {
            res = r;
            res.material = getMaterial(i, res.texCoords);
        }
    }
    return res;
}

// GPU random
uint seed = uint(
    uint((gl_FragCoord.x * 0.5 + 0.5) * screenWidth)  * uint(1973) + 
    uint((gl_FragCoord.y * 0.5 + 0.5) * screenHeight) * uint(9277) + 
    uint(counter) * uint(26699)) | uint(1);
uint hash(inout uint seed)
{
    seed = uint(seed ^ uint(61)) ^ uint(seed >> uint(16));
    seed *= uint(9);
    seed = seed ^ (seed >> 4);
    seed *= uint(0x27d4eb2d);
    seed = seed ^ (seed >> 15);
    return seed;
}
float rand()
{
    return float(hash(seed)) / 4294967296.0;
}
vec3 SampleHemisphere() {
    float z = rand();
    float r = max(0, sqrt(1.0 - z*z));
    float phi = 2.0 * PI * rand();
    return vec3(r * cos(phi), r * sin(phi), z);
}

// One ray tracing
vec3 tracing(HitResult hit, int maxBounce)
{
    vec3 Lo = vec3(0);
    vec3 nowL = vec3(1);

    for (int bounce = 0; bounce < maxBounce; bounce++) {
        vec3 wi = hit.TBN * SampleHemisphere();

        Ray randomRay;
        randomRay.startPoint = hit.hitPoint;
        randomRay.direction = wi;
        HitResult newHit = hitArray(randomRay, 0, nTriangles-1);

        float pdf = 1.0 / (2.0 * PI);
        float cosine_o = max(0, dot(-hit.viewDir, hit.normal));
        float cosine_i = max(0, dot(randomRay.direction, hit.normal));
        vec3 f_r = hit.material.diffuse / PI;

        if (!newHit.isHit) {
            break;
        }
        
        vec3 Le = newHit.material.emissive;
        Lo += nowL * Le * f_r * cosine_i / pdf;
        
        hit = newHit;
        nowL *= f_r * cosine_i / pdf;
    }
    return Lo;
}

void main()
{
    vec3 color = vec3(0.0);

    for (int i = 0; i < SAMPLENUM; i ++) {
        vec2 fragCoord = vec2(gl_FragCoord.x + rand() - 0.5, gl_FragCoord.y + rand() - 0.5);
        vec4 ndc = vec4(
            (fragCoord.x/screenWidth - 0.5) * 2.0,
            (fragCoord.y/screenHeight - 0.5) * 2.0,
            (gl_FragCoord.z - 0.5) * 2.0,
            1.0
        );
        vec4 clip = inverse(view) * inverse(projection) * ndc; 
        vec3 coord_pos = (clip / clip.w).xyz;

        Ray ray;
        ray.startPoint = viewPos;
        vec3 dir = coord_pos - ray.startPoint;
        ray.direction = normalize(dir);
        HitResult first = hitArray(ray, 0, nTriangles-1);

        if (first.isHit) {
            color += first.material.emissive + tracing(first, 2);
            // color += vec4(dot(first.normal, -first.viewDir), 0, 0, 1.0);
        }
        else {
            color += vec3(1.0);
        }
    }

    vec2 lastFrameTexCoord = vec2(gl_FragCoord.x/screenWidth, gl_FragCoord.y/screenHeight);
    FragColor = vec4(color / SAMPLENUM + texture2D(lastFrame, lastFrameTexCoord).rgb, 1.0);
    // FragColor = vec4(calcLerp(vec3(-1, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0), vec3(0.5, 0.5, 0)), 1.0);
    

    // Ray ray;
    // ray.startPoint = viewPos;
    // vec3 dir = coord_pos - ray.startPoint;
    // ray.direction = normalize(dir);
    // HitResult res = hitArray(ray, 0, nTriangles-1);

    // if (res.isHit)
    // // FragColor = vec4(gl_FragCoord.z / 4, 0.0, 0.0, 1.0);
    // FragColor = vec4(dot(res.normal, -res.viewDir), 0, 0, 1.0);
    // else FragColor = vec4(1.0, 1.0, 1.0, 1.0);
}