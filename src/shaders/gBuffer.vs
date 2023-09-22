#version 330 core

layout (location = 0) in vec3 a_Position;
layout (location = 1) in vec3 a_Normal;
layout (location = 2) in vec2 a_Texcoord;
layout (location = 3) in float a_MatID;

out vec3 v_Position;
out vec3 v_Normal;
out vec2 v_Texcoord;
out float v_MatID;
out vec4 nowVertex;
out vec4 prevVertex;

uniform mat4 prevPVMatrix;
uniform mat4 PVMatrix;

void main()
{
    prevVertex = prevPVMatrix * vec4(a_Position, 1.0f);
    nowVertex = PVMatrix * vec4(a_Position, 1.0f);

    gl_Position = PVMatrix * vec4(a_Position, 1.0f);

    v_Position = a_Position;
    v_Normal = a_Normal;
    v_Texcoord = a_Texcoord;
    v_MatID = a_MatID;
}