#version 330 core
layout (location = 0) out vec4 FragColor;

in vec3 FragCoord;

uniform sampler2D fragPos;
uniform sampler2D fragNorm;
uniform sampler2D fragAlbedo;
uniform sampler2D motionVecRoughness;

uniform sampler2D nowFrame;
uniform sampler2D lastFrame;

uniform bool isRealTime;
uniform int frameCounter;

#define SIGMA_TEX 3.0
#define SIGMA_POS 10.0
#define SIGMA_COLOR 5
#define SIGMA_ALBEDO 0.1
#define SIGMA_NORM 0.2

vec3 RGB2YCoCgR(vec3 rgbColor)
{
    vec3 YCoCgRColor;

    YCoCgRColor.y = rgbColor.r - rgbColor.b;
    float temp = rgbColor.b + YCoCgRColor.y / 2;
    YCoCgRColor.z = rgbColor.g - temp;
    YCoCgRColor.x = temp + YCoCgRColor.z / 2;

    return YCoCgRColor;
}
vec3 YCoCgR2RGB(vec3 YCoCgRColor)
{
    vec3 rgbColor;

    float temp = YCoCgRColor.x - YCoCgRColor.z / 2;
    rgbColor.g = YCoCgRColor.z + temp;
    rgbColor.b = temp - YCoCgRColor.y / 2;
    rgbColor.r = rgbColor.b + YCoCgRColor.y;

    return rgbColor;
}
float Luminance(vec3 color)
{
    return 0.25 * color.r + 0.5 * color.g + 0.25 * color.b;
}
vec3 ToneMap(vec3 color)
{
    return color / (1 + Luminance(color));
}
vec3 UnToneMap(vec3 color)
{
    return color / (1 - Luminance(color));
}

// clip color in aabbMin and aabbMax
vec3 clip(vec3 aabbMin, vec3 aabbMax, vec3 color)
{
    color = RGB2YCoCgR(ToneMap(color));
    vec3 p_clip = 0.5 * (aabbMax + aabbMin);
    vec3 e_clip = 0.5 * (aabbMax - aabbMin);

    vec3 v_clip = color - p_clip;
    vec3 v_unit = v_clip.xyz / e_clip;
    vec3 a_unit = abs(v_unit);
    float ma_unit = max(a_unit.x, max(a_unit.y, a_unit.z));

    if (ma_unit > 1.0) {
        color = p_clip + v_clip / ma_unit;
    }
    color = UnToneMap(YCoCgR2RGB(color));
    return color;
}

void main()
{
    vec2 uv = (FragCoord.xy + 1) / 2;

    // bilateral filter
    int R = 31;
    float sumW = 0;
    vec3 centerColor = texture(nowFrame, uv).rgb;
    vec2 texelSize = 1.0 / textureSize(nowFrame, 0).xy;
    vec3 pos = texture(fragPos, uv).rgb;
    vec3 norm = texture(fragNorm, uv).rgb;
    vec3 albedo = texture(fragAlbedo, uv).rgb;
    vec3 m1 = vec3(0);
    vec3 m2 = vec3(0);
    vec3 color = vec3(0);
    for (int i = -R; i <= R; i += 2) {
        for (int j = -R; j <= R; j += 2) {
            vec2 newUV = uv + vec2(i, j) * texelSize;
            vec3 newPos = texture(fragPos, newUV).rgb;
            vec3 newNorm = texture(fragNorm, newUV).rgb;
            vec3 newAlbedo = texture(fragAlbedo, newUV).rgb;
            vec3 newColor = texture(nowFrame, newUV).rgb;
            float smoothness = 1 - texture(motionVecRoughness, newUV).b;
            float wexp = 0;
            wexp -= pow(i * i + j * j, smoothness) / (2 * SIGMA_TEX * SIGMA_TEX);
            wexp -= pow(length(newPos - pos), 2) / (2 * SIGMA_POS * SIGMA_POS);
            wexp -= pow(length(newColor - centerColor), 2) / (2 * SIGMA_COLOR * SIGMA_COLOR);
            wexp -= pow(length(newAlbedo - albedo), 2) / (2 * SIGMA_ALBEDO * SIGMA_ALBEDO);
            wexp -= pow(length(newNorm - norm), 2) / (2 * SIGMA_NORM * SIGMA_NORM);
            float w = exp(wexp);
            sumW += w;
            color += w * newColor;

            newColor = RGB2YCoCgR(ToneMap(newColor));
            m1 += newColor;
            m2 += newColor * newColor;
        }
    }
    color /= sumW;
    // outlier removal
    int N = (2 * R + 1) * (2 * R + 1);
    float VarianceClipGamma = 1.0f;
    vec3 mu = m1 / N;
    vec3 sigma = sqrt(abs(m2 / N - mu * mu));
    vec3 aabbMin = mu - VarianceClipGamma * sigma;
    vec3 aabbMax = mu + VarianceClipGamma * sigma;
    color = clip(aabbMin, aabbMax, color);

    // temporal
    vec2 mVec = texture(motionVecRoughness, uv).rg;
    vec2 prevUV = uv - mVec;
    vec3 lastColor = texture(lastFrame, prevUV).rgb;
    if (prevUV.x < 0 || prevUV.x > 1 || prevUV.y < 0 || prevUV.y > 1) {
        lastColor = color;
    }
    else {
        lastColor = clip(aabbMin, aabbMax, lastColor);
    }

    if (isRealTime) {
        color = mix(lastColor, color, 0.1);
    }
    else {
        color = length(color) >= 0.0 ? mix(lastColor, color, 1.0/frameCounter) : lastColor;
    }

    FragColor = vec4(color, 1.0);
}