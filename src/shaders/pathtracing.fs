#version 330 core
layout (location = 0) out vec4 FragColor;

in vec3 FragCoord;
in mat4 inverse_projection;
in mat4 inverse_view;

uniform bool isRealTime;
uniform float screenWidth;
uniform float screenHeight;

uniform int frameCounter;

uniform sampler2DArray texturesArray;
uniform samplerBuffer triangles;
uniform int nTriangles;
uniform samplerBuffer bvhNodes;
uniform int nNodes;

uniform sampler2D lastFrame;
uniform sampler2D fragPos;
uniform sampler2D hdrMap;
uniform vec3 viewPos;

#define STACKSIZE 64
#define SAMPLENUM 1
#define PI 3.14159265359

struct Material {
    vec3 ambient;       // texture: (-id-1, width, height) color: (r, g, b)
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

struct Triangle {
    vec3 p[3];
    vec3 n[3];
    vec2 texCoords[3];
    mat3 TBN[3];
    vec3 normal;        // texture: (-id-1, width, height) none: (-1, -1, -1)
    vec3 height;        // texture: (-id-1, width, height) none: (-1, -1, -1)
};

struct BVHNode {
    int left, right;
    int n;
    int index1, index2;
    int parent;
    vec3 AA, BB;
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

// GPU random
uint seed = uint(
    uint(gl_FragCoord.x) * uint(1973) + 
    uint(gl_FragCoord.y) * uint(9277) + 
    uint(frameCounter) * uint(26699)) | uint(1);
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
float sqr(float x) { 
    return x*x; 
}
vec3 SampleHemisphere() {
    float z = rand();
    float r = max(0, sqrt(1.0 - z*z));
    float phi = 2.0 * PI * rand();
    return vec3(r * cos(phi), r * sin(phi), z);
}
vec3 SampleGTR2(vec3 V, float alpha, mat3 TBN)
{
    float xi_1 = rand();
    float xi_2 = rand();
    float phi_h = 2.0 * PI * xi_1;
    float sin_phi_h = sin(phi_h);
    float cos_phi_h = cos(phi_h);

    float cos_theta_h = sqrt((1.0-xi_2)/(1.0+(alpha*alpha-1.0)*xi_2));
    float sin_theta_h = sqrt(max(0.0, 1.0 - cos_theta_h * cos_theta_h));

    vec3 H = vec3(sin_theta_h*cos_phi_h, sin_theta_h*sin_phi_h, cos_theta_h);
    H = TBN * H;

    vec3 L = reflect(-V, H);

    return L;
}

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

vec3 getColor(int pos, vec2 texCoords)
{
    vec3 t = texelFetch(triangles, pos).xyz;
    if (t.x < 0) {
        return texture2DArray(texturesArray, vec3(texCoords * t.yz, - t.x - 1)).rgb;
    }
    else {
        return t;
    }
}

Material getMaterial(int i, vec2 texCoords)
{
    int offset = i * 29;
    Material m;
    vec3 tmp;

    m.ambient = getColor(offset + 19, texCoords);
    m.diffuse = getColor(offset + 20, texCoords);
    m.specular = getColor(offset + 21, texCoords);
    m.emissive = getColor(offset + 22, texCoords);
    m.shininess = getColor(offset + 23, texCoords).x;
    m.roughness = 1 - m.shininess / 1000;
    m.metallic = getColor(offset + 24, texCoords).x;
    m.refracti = getColor(offset + 25, texCoords).x;
    m.opacity = getColor(offset + 26, texCoords).x;
    m.transmission = getColor(offset + 27, texCoords);
    m.anisotropy = getColor(offset + 28, texCoords).x;

    return m;
}

BVHNode getBVHNode(int i)
{
    int offset = i * 4;
    BVHNode node;

    node.left = int(texelFetch(bvhNodes, offset).x);
    node.right = int(texelFetch(bvhNodes, offset).y);
    node.n = int(texelFetch(bvhNodes, offset).z);
    node.index1 = int(texelFetch(bvhNodes, offset + 1).x);
    node.index2 = int(texelFetch(bvhNodes, offset + 1).y);
    node.parent = int(texelFetch(bvhNodes, offset + 1).z);
    node.AA = texelFetch(bvhNodes, offset + 2).xyz;
    node.BB = texelFetch(bvhNodes, offset + 3).xyz;

    return node;
}

HitResult hitTriangle(Triangle triangle, Ray ray)
{
    HitResult res;
    res.distance = 0;
    res.isHit = false;
    res.isInside = false;

    vec3 O = ray.startPoint;
    vec3 D = ray.direction;
    vec3 E1 = triangle.p[1] - triangle.p[0];
    vec3 E2 = triangle.p[2] - triangle.p[0];
    vec3 S1 = cross(D, E2);

    float SE = dot(S1, E1);
    // parallel
    if (abs(SE) < 0.000001f) return res;

    float invdet = 1 / SE;
    vec3 S = O - triangle.p[0];

    float b = dot(S1, S) * invdet;
    // out of triangle
    if (b < 0 || b > 1) return res;

    vec3 S2 = cross(S, E1);
    float c = dot(S2, D) / SE;
    // out of triangle
    if (c < 0 || b + c > 1) return res;

    float t = dot(S2, E2) / SE;
    if (t < 0.0005f) return res;

    float a = 1 - b - c;
    res.isHit = true;
    res.distance = t;
    res.hitPoint = O + t * D;
    res.viewDir = D;
    res.texCoords = a * triangle.texCoords[0] + b * triangle.texCoords[1] + c * triangle.texCoords[2];
    res.TBN = a * triangle.TBN[0] + b * triangle.TBN[1] + c * triangle.TBN[2];

    // count normal and TBN
    if (triangle.normal.y < 0 && triangle.normal.z < 0) {
        res.normal = a * triangle.n[0] + b * triangle.n[1] + c * triangle.n[2];
    }
    else {
        // normal Map
        res.normal = texture2DArray(texturesArray, vec3(res.texCoords * triangle.normal.yz, triangle.normal.x)).rgb;
        res.normal = res.normal * 2.0 - 1.0;
        res.normal = res.TBN * res.normal;
    }
    
    if (dot(res.normal, D) > 0) res.normal = -res.normal;
    res.normal = normalize(res.normal);

    vec3 tangent = normalize(cross(res.TBN * vec3(0, 1, 0), res.normal));
    vec3 bitangent = normalize(cross(res.normal, tangent));
    res.TBN = mat3(tangent, bitangent, res.normal);

    return res;
}

HitResult hitArray(Ray ray, int l, int r)
{
    HitResult res;
    res.isHit = false;
    for (int i = l; i <= r; i ++) {
        Triangle triangle = getTriangle(i);
        HitResult newhit = hitTriangle(triangle, ray);
        if (newhit.isHit && (!res.isHit || newhit.distance < res.distance)) {
            res = newhit;
            res.material = getMaterial(i, res.texCoords);
        }
    }
    return res;
}

// if hit AABB, return distance
float hitAABB(Ray r, vec3 AA, vec3 BB)
{
    vec3 invdir = 1.0 / r.direction;

    vec3 f = (BB - r.startPoint) * invdir;
    vec3 n = (AA - r.startPoint) * invdir;

    vec3 tmax = max(f, n);
    vec3 tmin = min(f, n);

    float t1 = min(tmax.x, min(tmax.y, tmax.z));
    float t0 = max(tmin.x, max(tmin.y, tmin.z));

    return (t1 >= t0) ? ((t0 > 0.0) ? (t0) : (t1)) : (-1);
}

// search in BVH tree
HitResult hitBVH(Ray ray)
{
    HitResult res;
    res.isHit = false;

    int stack[STACKSIZE];
    int sp = 0;

    stack[sp++] = 0;
    while (sp > 0) {
        int top = stack[--sp];
        BVHNode node = getBVHNode(top);
        
        if (node.n > 0) {
            HitResult r = hitArray(ray, node.index1, node.index2);
            if (r.isHit && (!res.isHit || r.distance < res.distance)) res = r;
            continue;
        }
        
        float d1, d2;
        if (node.left > 0) {
            BVHNode leftNode = getBVHNode(node.left);
            d1 = hitAABB(ray, leftNode.AA, leftNode.BB);
        }
        if (node.right > 0) {
            BVHNode rightNode = getBVHNode(node.right);
            d2 = hitAABB(ray, rightNode.AA, rightNode.BB);
        }

        if (d1 > 0 && d2 > 0) {
            if (d1 < d2) {
                stack[sp++] = node.right;
                stack[sp++] = node.left;
            }
            else {
                stack[sp++] = node.left;
                stack[sp++] = node.right;
            }
        }
        else if (d1 > 0) {
            stack[sp++] = node.left;
        }
        else if (d2 > 0) {
            stack[sp++] = node.right;
        }
    }

    return res;
}

// sample hdr
vec2 sampleSphericalMap(vec3 v) {
    vec2 uv = vec2(atan(v.z, v.x), asin(v.y));
    uv /= vec2(2.0 * PI, PI);
    uv += 0.5;
    uv.y = 1.0 - uv.y;
    return uv;
}
vec3 sampleHdr(vec3 v) {
    vec2 uv = sampleSphericalMap(normalize(v));
    vec3 color = texture2D(hdrMap, uv).rgb;
    //color = min(color, vec3(10));
    return color;
}

// count BRDF
float SchlickFresnel(float u) {
    float m = clamp(1-u, 0, 1);
    float m2 = m*m;
    return m2 * m2 * m;
}
float GTR2_aniso(float NdotH, float HdotX, float HdotY, float ax, float ay) {
    return 1 / (PI * ax*ay * sqr( sqr(HdotX/ax) + sqr(HdotY/ay) + NdotH*NdotH ));
}
float smithG_GGX_aniso(float NdotV, float VdotX, float VdotY, float ax, float ay) {
    return 1 / (NdotV + sqrt( sqr(VdotX*ax) + sqr(VdotY*ay) + sqr(NdotV) ));
}
vec3 BRDF(vec3 V, vec3 N, vec3 L, mat3 TBN, Material mat)
{
    float NdotL = dot(N, L);
    float NdotV = dot(N, V);
    if(NdotL < 0 || NdotV < 0) return vec3(0);

    vec3 H = normalize(V + L);
    float NdotH = dot(N, H);
    float LdotH = dot(L, H);
    vec3 X = TBN * vec3(1, 0, 0);
    vec3 Y = TBN * vec3(0, 1, 0);

    // diffuse
    float Fd90 = 0.5 + 2.0 * LdotH * LdotH * mat.roughness;
	float FL = SchlickFresnel(NdotL);
	float FV = SchlickFresnel(NdotV);
	float Fd = mix(1.0, Fd90, FL) * mix(1.0, Fd90, FV);
	vec3 diffuse = Fd * mat.diffuse / PI;

    // specular
    float aspect = sqrt(1.0 - mat.anisotropy * 0.9);
    float ax = max(0.001, sqr(mat.roughness)/aspect);
    float ay = max(0.001, sqr(mat.roughness)*aspect);
    float F0 = (mat.refracti - 1) * (mat.refracti - 1) / (mat.refracti + 1) / (mat.refracti + 1);
    vec3 spec = mix(F0 * mat.specular, mat.diffuse, mat.metallic);

    float Ds = GTR2_aniso(NdotH, dot(H, X), dot(H, Y), ax, ay);
    float Gs = smithG_GGX_aniso(NdotV, dot(V, X), dot(V, Y), ax, ay);
    Gs *= smithG_GGX_aniso(NdotL, dot(L, X), dot(L, Y), ax, ay);
    vec3 Fs = mix(spec, vec3(1), SchlickFresnel(LdotH));
    vec3 specular = Ds * Gs * Fs;

	return (1.0 - mat.metallic) * diffuse + specular;
}
float BRDF_Pdf(vec3 V, vec3 N, vec3 L, mat3 TBN, Material material) {
    float NdotL = dot(N, L);
    float NdotV = dot(N, V);
    if (NdotL < 0 || NdotV < 0) return 0;

    vec3 H = normalize(L + V);
    float NdotH = dot(N, H);
    float LdotH = dot(L, H);
    vec3 X = TBN * vec3(1, 0, 0);
    vec3 Y = TBN * vec3(0, 1, 0);
     
    float aspect = sqrt(1.0 - material.anisotropy * 0.9);
    float ax = max(0.001, sqr(material.roughness)/aspect);
    float ay = max(0.001, sqr(material.roughness)*aspect);
    float Ds = GTR2_aniso(NdotH, dot(H, X), dot(H, Y), ax, ay);

    float pdf_diffuse = NdotL / PI;
    float pdf_specular = Ds * NdotH / (4.0 * dot(L, H));

    float r_diffuse = (1.0 - material.metallic);
    float r_specular = 1.0;
    float r_sum = r_diffuse + r_specular;

    float p_diffuse = r_diffuse / r_sum;
    float p_specular = r_specular / r_sum;

    float pdf = p_diffuse * pdf_diffuse + p_specular * pdf_specular;

    pdf = max(1e-10, pdf);
    return pdf;
}
vec3 SampleL(vec3 V, vec3 N, mat3 TBN, Material material)
{
    float alpha = max(0.001, sqr(material.roughness));
    
    float r_diffuse = (1.0 - material.metallic);
    float r_specular = 1.0;
    float r_sum = r_diffuse + r_specular;

    float p_diffuse = r_diffuse / r_sum;
    float p_specular = r_specular / r_sum;

    float p = rand();

    if(p <= p_diffuse) {
        return TBN * SampleHemisphere();
    } 
    else if(p <= p_diffuse + p_specular) {    
        return SampleGTR2(V, alpha, TBN);
    }
    return vec3(0, 1, 0);
}

// One ray tracing
vec3 tracing(HitResult hit, int maxBounce)
{
    vec3 Lo = vec3(0);
    vec3 nowL = vec3(1);

    for(int bounce=0; bounce<maxBounce; bounce++) {
        vec3 V = -hit.viewDir;
        vec3 N = hit.normal;

        vec3 L = SampleL(V, N, hit.TBN, hit.material); 
        float NdotL = dot(N, L);
        if (NdotL <= 0.0) break;

        Ray randomRay;
        randomRay.startPoint = hit.hitPoint;
        randomRay.direction = L;
        HitResult newHit = hitBVH(randomRay);

        vec3 f_r = BRDF(V, N, L, hit.TBN, hit.material);
        float pdf = BRDF_Pdf(V, N, L, hit.TBN, hit.material);
        if (pdf <= 0.0) break;

        if(!newHit.isHit) {
            vec3 skyColor = sampleHdr(randomRay.direction);
            Lo += nowL * skyColor * f_r * NdotL / pdf;
            break;
        }
        
        vec3 Le = newHit.material.emissive;
        Lo += nowL * Le * f_r * NdotL / pdf;             

        hit = newHit;
        nowL *= f_r * NdotL / pdf;
    }
    
    return Lo;
}

void main()
{
    vec3 color = vec3(0.0);

    for (int i = 0; i < SAMPLENUM; i ++) {
        vec2 randxy = vec2((rand() - 0.5) / screenWidth, (rand() - 0.5) / screenHeight);
        vec4 ndc = vec4(FragCoord.xy + randxy, FragCoord.z, 1.0);
        vec4 clip = inverse_view * inverse_projection * ndc; 
        vec3 coord_pos = (clip / clip.w).xyz;
        
        Ray ray;
        ray.startPoint = viewPos;
        vec3 dir = coord_pos - ray.startPoint;
        ray.direction = normalize(dir);
        HitResult first = hitBVH(ray);

        if (first.isHit) {
            color += first.material.emissive + tracing(first, 2);
        }
        else {
            color += sampleHdr(ray.direction);
        }
    }

    color /= SAMPLENUM;
    vec3 lastColor = texture(lastFrame, (FragCoord.xy + 1) / 2).rgb;
    if (!isRealTime) {
        color = length(color) > 0.0 ? mix(lastColor, color, 1/log(9+frameCounter)) : lastColor;
    }
    FragColor = vec4(color, 1.0);
}