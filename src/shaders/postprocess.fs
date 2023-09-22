#version 330 core
layout (location = 0) out vec4 FragColor;

in vec3 FragCoord;

uniform sampler2D frame;

vec3 toneMapping(in vec3 c, float limit) {
    float luminance = 0.3*c.x + 0.6*c.y + 0.1*c.z;
    return c * 1.0 / (1.0 + luminance / limit);
}

void main()
{
    vec2 uv = (FragCoord.xy + 1) / 2;
    vec2 texelSize = 1.0 / textureSize(frame, 0).xy;

    vec3 color = texture(frame, uv).rgb;
    color = toneMapping(color, 1.5);
    color = pow(color, vec3(1.0 / 2.2));

    // int R = 3;
    // vec3 m1 = vec3(0);
    // vec3 m2 = vec3(0);
    // for (int i = -R; i <= R; i ++) {
    //     for (int j = -R; j <= R; j ++) {
    //         vec2 newUV = uv + vec2(i, j) * texelSize;
    //         vec3 newColor = texture(frame, newUV).rgb;
    //         newColor = toneMapping(newColor, 1.5);
    //         newColor = pow(newColor, vec3(1.0 / 2.2));
    //         m1 += newColor;
    //         m2 += newColor * newColor;
    //     }
    // }

    // int N = (2 * R + 1) * (2 * R + 1);
    // float VarianceClipGamma = 1.0f;
    // vec3 mu = m1 / N;
    // vec3 sigma = sqrt(abs(m2 / N - mu * mu));
    // vec3 aabbMin = mu - VarianceClipGamma * sigma;
    // vec3 aabbMax = mu + VarianceClipGamma * sigma;

    // if (color.x < aabbMin.x || color.y < aabbMin.y || color.z < aabbMin.z) {
    //     color = aabbMin;
    // }
    // if (color.x > aabbMax.x || color.y > aabbMax.y || color.z > aabbMax.z) {
    //     color = aabbMax;
    // }

    FragColor = vec4(color, 1.0);
}