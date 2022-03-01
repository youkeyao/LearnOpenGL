#version 330 core
layout (location = 0) out vec3 gPosition;

in vec3 FragPos;

void main()
{    
    gPosition = FragPos;
}