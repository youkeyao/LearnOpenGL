#version 330 core
layout (location = 0) out vec4 FragColor;

in vec3 FragCoord;

uniform sampler2D fragPos;
uniform sampler2D nowFrame;

uniform int level = 0;
uniform int support = 2;
uniform float sigmaP = 1.0;
uniform float sigmaN = 128.0;
uniform float sigmaL = 4.0;

const float epsilon = 0.00001;
const float h[25] = {1.0/256.0, 1.0/64.0, 3.0/128.0, 1.0/64.0, 1.0/256.0,
                              1.0/64.0,  1.0/16.0, 3.0/32.0,  1.0/16.0, 1.0/64.0,
                              3.0/128.0, 3.0/32.0, 9.0/64.0,  3.0/32.0, 3.0/128.0,
                              1.0/64.0,  1.0/16.0, 3.0/32.0,  1.0/16.0, 1.0/64.0,
                              1.0/256.0, 1.0/64.0, 3.0/128.0, 1.0/64.0, 1.0/256.0};
const float gaussKernel[9] = {1.0/16.0, 1.0/8.0, 1.0/16.0, 1.0/8.0, 1.0/4.0, 1.0/8.0, 1.0/16.0, 1.0/8.0, 1.0/16.0};

float luma(vec3 c){
    return dot(c, vec3(0.2126, 0.7152, 0.0722));
}

void main()
{
    // vec2 uv = (FragCoord.xy + 1) / 2;

    // vec3 pPosition = texture(fragPos, uv).rgb;

    // vec3 pColor = texture(nowFrame, uv).rgb;
    // float pLuminance = luma(pColor);

    // vec2 texelSize = 1.0 / textureSize(nowFrame, 0).xy;
    // int step = 1 << level;

    // vec3 c = vec3(0.0);
    // float v = 0.0;
    // float weights = 0.0;

    // for (int offsetx = -support; offsetx <= support; offsetx++) {
    //     for (int offsety = -support; offsety <= support; offsety++) {
    //         vec2 loc = uv + vec2(step * offsetx * texelSize.x, step * offsety * texelSize.y);

    //         // if (eps_equal(pMeshID, qMeshID)) {
    //             vec3 qPosition = texture(fragPos, loc).rgb;

    //             vec3 qColor = texture(nowFrame, loc).rgb;
    //             float qVariance = texture(nowFrame, loc).a;
    //             float qLuminance = luma(qColor);

    //             vec3 t = pPosition - qPosition;
    //             float dist2 = dot(t, t) + t.z * t.z;
    //             float wp = min(exp(-(dist2)/sigmaP), 1.0);

    //             float wn = sigmaN;

    //             float gvl = 0.001;
    //             for (int y0 = -1; y0 <= 1; y0++) {
    //                 for (int x0 = -1; x0 <= 1; x0++) {
    //                     gvl += gaussKernel[x0 + 3*y0 + 4] * texture(nowFrame, loc + vec2(x0, y0) * texelSize).a;
    //                 }
    //             }
    //             float wl = min(1.0, exp(-abs(pLuminance - qLuminance) / (sigmaL * sqrt(gvl) + epsilon)));

    //             float w = wp * wn * wl;
    //             float weight = h[5*(offsety + support) + offsetx + support] * w;

    //             c += weight * qColor;
    //             v += weight * weight * qVariance;
    //             weights += weight;
    //         // }
    //     }
    // }

    // if (weights > epsilon) {
    //     FragColor.rgb = c / weights;
    //     FragColor.a = v / (weights * weights);
    // } else {
    //     FragColor = texture(nowFrame, uv).rgba;
    // }

    FragColor = vec4(texture(nowFrame, (FragCoord.xy + 1) / 2).rgb, 1.0);
}