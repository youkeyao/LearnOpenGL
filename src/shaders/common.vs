#version 330 core
layout (location = 0) in vec3 position;

out vec3 FragCoord;
out mat4 inverse_projection;
out mat4 inverse_view;

uniform mat4 projection;
uniform mat4 view;

void main()
{
    gl_Position = vec4(position, 1.0);
    FragCoord = position;
    inverse_projection = inverse(projection);
    inverse_view = inverse(view);
}