#ifndef RENDERPASS_H
#define RENDERPASS_H

#include <glad/glad.h>

#include <vector>
using namespace std;

class RenderPass
{
public:
    unsigned int VAO, VBO, FBO = 0;
    vector<unsigned int> attachments;
    RenderPass(bool finalPass)
    {
        if (!finalPass) glGenFramebuffers(1, &FBO);
        GLfloat quadVertices[] = {
            // Positions        // Texture Coords
            -1.0f, 1.0f, 0.0f, 0.0f, 1.0f,
            -1.0f, -1.0f, 0.0f, 0.0f, 0.0f,
            1.0f, 1.0f, 0.0f, 1.0f, 1.0f,
            1.0f, -1.0f, 0.0f, 1.0f, 0.0f,
        };
        // Setup plane VAO
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), (GLvoid*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
        glBindVertexArray(0);
    }

    void bindData(unsigned int type, int width, int height)
    {
        glBindFramebuffer(GL_FRAMEBUFFER, FBO);
        if (type == GL_DEPTH_ATTACHMENT) {
            unsigned int rboDepth;
            glGenRenderbuffers(1, &rboDepth);
            glBindRenderbuffer(GL_RENDERBUFFER, rboDepth);
            glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, width, height);
            glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rboDepth);
        }
        else {
            unsigned int colorAttachment;
            glGenTextures(1, &colorAttachment);
            glActiveTexture(GL_TEXTURE0 + colorAttachment);
            glBindTexture(GL_TEXTURE_2D, colorAttachment);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, width, height, 0, GL_RGB, GL_FLOAT, NULL);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glFramebufferTexture2D(GL_FRAMEBUFFER, type, GL_TEXTURE_2D, colorAttachment, 0);
            glActiveTexture(0);
            attachments.push_back(colorAttachment);
        }
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    void activeFramebuffer(vector<unsigned int> texPass = {})
    {
        glBindFramebuffer(GL_FRAMEBUFFER, FBO);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        for (int i = 0; i < texPass.size(); i ++) {
            glActiveTexture(GL_TEXTURE0 + texPass[i]);
            glBindTexture(GL_TEXTURE_2D, texPass[i]);
        }
    }

    void draw()
    {
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glBindVertexArray(0);
    }
};

#endif