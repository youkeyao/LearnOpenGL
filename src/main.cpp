#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <utility/Shader.h>
#include <utility/Camera.h>
#include <utility/Scene.h>
#include <utility/RenderPass.h>

#include <iostream>

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow *window);

// settings
const unsigned int SCR_WIDTH = 600;
const unsigned int SCR_HEIGHT = 600;
bool isRealTime = true;
bool isRealTimeRelease = true;
bool isFocus = true;
bool isFocusRelease = true;

// camera
Camera camera(glm::vec3(0.0f, -4.0f, 7.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

// timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;
int frameCounter = 0;

// dir
std::string shadersDir = "../src/shaders/";
std::string resourcesDir = "../resources/";

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);

    // build and compile shaders
    // -------------------------
    Shader shaderGeometryPass((shadersDir + "gBuffer.vs").c_str(), (shadersDir + "gBuffer.fs").c_str());
    Shader shaderPathTracingPass((shadersDir + "common.vs").c_str(), (shadersDir + "pathtracing.fs").c_str());
    Shader shaderFilterPass((shadersDir + "common.vs").c_str(), (shadersDir + "filter.fs").c_str());
    Shader shaderOutputPass((shadersDir + "common.vs").c_str(), (shadersDir + "output.fs").c_str());
    Shader shaderPostProcessPass((shadersDir + "common.vs").c_str(), (shadersDir + "postprocess.fs").c_str());
    shaderPathTracingPass.use();
    shaderPathTracingPass.setFloat("screenWidth", SCR_WIDTH);
    shaderPathTracingPass.setFloat("screenHeight", SCR_HEIGHT);

    // load scene
    // -----------
    Scene scene(resourcesDir + "objects/testScene/testScene.obj", aiProcess_Triangulate | aiProcess_GenSmoothNormals | aiProcess_CalcTangentSpace | aiProcess_FlipUVs);
    scene.loadShader(shaderPathTracingPass);
    scene.loadShader(shaderGeometryPass);
    std::cout << "scene loaded" << std::endl;

    // load the HDR environment map
    // ---------------------------------
    stbi_set_flip_vertically_on_load(false);
    int width, height, nrComponents;
    float *data = stbi_loadf((resourcesDir + "textures/hdr/newport_loft.hdr").c_str(), &width, &height, &nrComponents, 0);
    unsigned int hdrTexture;
    if (data) {
        glGenTextures(1, &hdrTexture);
        glActiveTexture(GL_TEXTURE0 + hdrTexture);
        glBindTexture(GL_TEXTURE_2D, hdrTexture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, width, height, 0, GL_RGB, GL_FLOAT, data);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        glActiveTexture(0);

        stbi_image_free(data);
    }
    else
    {
        std::cout << "Failed to load HDR image." << std::endl;
    }

    // Init render pass
    // ---------------------------------
    RenderPass geometryPass(false);
    geometryPass.bindData(GL_DEPTH_ATTACHMENT, SCR_WIDTH, SCR_HEIGHT);
    geometryPass.bindData(GL_COLOR_ATTACHMENT0, SCR_WIDTH, SCR_HEIGHT);
    geometryPass.bindData(GL_COLOR_ATTACHMENT1, SCR_WIDTH, SCR_HEIGHT);
    geometryPass.bindData(GL_COLOR_ATTACHMENT2, SCR_WIDTH, SCR_HEIGHT);
    geometryPass.bindData(GL_COLOR_ATTACHMENT3, SCR_WIDTH, SCR_HEIGHT);
    RenderPass pathTracingPass(false);
    pathTracingPass.bindData(GL_COLOR_ATTACHMENT0, SCR_WIDTH, SCR_HEIGHT);
    RenderPass filterPass(false);
    filterPass.bindData(GL_COLOR_ATTACHMENT0, SCR_WIDTH, SCR_HEIGHT);
    RenderPass outputPass(false);
    outputPass.bindData(GL_COLOR_ATTACHMENT0, SCR_WIDTH, SCR_HEIGHT);
    RenderPass postProcessPass(true);

    glm::mat4 pvMatrix = glm::mat4(1.0f);

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        // --------------------
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        frameCounter ++;

        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glm::mat4 projection = glm::perspective(camera.Zoom, (GLfloat)SCR_WIDTH / (GLfloat)SCR_HEIGHT, 0.1f, 100.0f);
        glm::mat4 view = camera.GetViewMatrix();

        // 1. Geometry Pass: render scene's geometry/color data into gbuffer
        // ---------------------------------
        geometryPass.activeFramebuffer();
        shaderGeometryPass.use();
        shaderGeometryPass.setMat4("prevPVMatrix", pvMatrix);
        pvMatrix = projection * view;
        shaderGeometryPass.setMat4("PVMatrix", pvMatrix);
        scene.draw();
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        // 2. PathTracing Pass: calculate lighting by iterating over a screen filled quad pixel-by-pixel using the gbuffer's content.
        // ---------------------------------
        pathTracingPass.activeFramebuffer();
        shaderPathTracingPass.use();
        shaderPathTracingPass.setMat4("projection", projection);
        shaderPathTracingPass.setMat4("view", view);
        shaderPathTracingPass.setVec3("viewPos", camera.Position);
        shaderPathTracingPass.setInt("hdrMap", hdrTexture);
        shaderPathTracingPass.setInt("frameCounter", frameCounter);
        pathTracingPass.draw();
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        // 3. Filter Pass: denoise
        // ---------------------------------
        filterPass.activeFramebuffer(geometryPass.attachments);
        filterPass.activeFramebuffer(pathTracingPass.attachments);
        filterPass.activeFramebuffer(outputPass.attachments);
        shaderFilterPass.use();
        shaderFilterPass.setInt("fragPos", geometryPass.attachments[0]);
        shaderFilterPass.setInt("fragNorm", geometryPass.attachments[1]);
        shaderFilterPass.setInt("fragAlbedo", geometryPass.attachments[2]);
        shaderFilterPass.setInt("motionVecRoughness", geometryPass.attachments[3]);
        shaderFilterPass.setInt("nowFrame", pathTracingPass.attachments[0]);
        shaderFilterPass.setInt("lastFrame", outputPass.attachments[0]);
        shaderFilterPass.setInt("frameCounter", frameCounter);
        shaderFilterPass.setBool("isRealTime", isRealTime);
        filterPass.draw();
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        // 4. output Pass: output last frame
        // ---------------------------------
        outputPass.activeFramebuffer(filterPass.attachments);
        shaderOutputPass.use();
        shaderOutputPass.setInt("nowFrame", filterPass.attachments[0]);
        outputPass.draw();
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        // 5. Postprocess Pass: postprocess and render to screen
        // ---------------------------------
        postProcessPass.activeFramebuffer(outputPass.attachments);
        shaderPostProcessPass.use();
        shaderPostProcessPass.setInt("frame", outputPass.attachments[0]);
        postProcessPass.draw();
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        if (isFocusRelease) {
            isFocusRelease = false;
            isFocus = !isFocus;
            if (isFocus) {
                glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
                glfwSetCursorPosCallback(window, mouse_callback);
            }
            else {
                glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
                glfwSetCursorPosCallback(window, NULL);
            }
            std::cout << "isFocus: " << isFocus << std::endl;
        }
    }
    else isFocusRelease = true;

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
        if (isRealTimeRelease) {
            isRealTimeRelease = false;
            isRealTime = !isRealTime;
            frameCounter = 1;
            std::cout << "isRealTime: " << isRealTime << std::endl;
        }
    }
    else isRealTimeRelease = true;
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(static_cast<float>(yoffset));
}