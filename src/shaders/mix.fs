#version 330 core
layout (location = 0) out vec4 FragColor;

in vec3 FragCoord;

uniform sampler2D nowFrame;

void main()
{    
    FragColor = vec4(texture(nowFrame, (FragCoord.xy + 1) / 2).rgb, 1.0);
}