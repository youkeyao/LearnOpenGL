#ifndef SCENE_H
#define SCENE_H

#include <glad/glad.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include <stb_image.h>
#include <utility/Shader.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>

struct Texture {
    glm::vec3 id;
    std::string name;
};

struct Material {
    glm::vec3 ambient;          // texture: (-id-2, width, height) color: (r, g, b)
    glm::vec3 diffuse;
    glm::vec3 specular;
    glm::vec3 emissive;
    glm::vec3 shininess;        // texture: (-id-2, width, height) color: (shininess, 0, 0)
    glm::vec3 metallic;
    glm::vec3 refracti;
    glm::vec3 opacity;
    glm::vec3 transmission;
    glm::vec3 anisotropy;
};

struct Triangle {
    float v[3];             // index of Vertex array
};

struct Vertex {
    glm::vec3 pos;
    glm::vec3 normal;           // texture: (-id-2, width, height) normal: (n1, n2, n3)
    glm::vec2 texCoords;
    float materialID;
    glm::vec3 tangent;
    glm::vec3 bitangent;
};

struct BVHNode {
    float left, right;          // index of children
    float n;                    // number of triangles in node
    float index1, index2;       // left and right index of triangles
    float parent;
    glm::vec3 AA, BB;
};

// help function
// ------------------
bool cmpx(const std::vector<Vertex>& vertices, const Triangle& t1, const Triangle& t2)
{
    glm::vec3 center1 = (vertices[t1.v[0]].pos + vertices[t1.v[1]].pos + vertices[t1.v[2]].pos);
    glm::vec3 center2 = (vertices[t2.v[0]].pos + vertices[t2.v[1]].pos + vertices[t2.v[2]].pos);
    return center1.x < center2.x;
}
bool cmpy(const std::vector<Vertex>& vertices, const Triangle& t1, const Triangle& t2)
{
    glm::vec3 center1 = (vertices[t1.v[0]].pos + vertices[t1.v[1]].pos + vertices[t1.v[2]].pos);
    glm::vec3 center2 = (vertices[t2.v[0]].pos + vertices[t2.v[1]].pos + vertices[t2.v[2]].pos);
    return center1.y < center2.y;
}
bool cmpz(const std::vector<Vertex>& vertices, const Triangle& t1, const Triangle& t2)
{
    glm::vec3 center1 = (vertices[t1.v[0]].pos + vertices[t1.v[1]].pos + vertices[t1.v[2]].pos);
    glm::vec3 center2 = (vertices[t2.v[0]].pos + vertices[t2.v[1]].pos + vertices[t2.v[2]].pos);
    return center1.z < center2.z;
}
glm::vec3 minVec3(glm::vec3 v1, glm::vec3 v2, glm::vec3 v3)
{
    float minx = std::min(v1.x, std::min(v2.x, v3.x));
    float miny = std::min(v1.y, std::min(v2.y, v3.y));
    float minz = std::min(v1.z, std::min(v2.z, v3.z));
    return glm::vec3(minx, miny, minz);
}
glm::vec3 maxVec3(glm::vec3 v1, glm::vec3 v2, glm::vec3 v3)
{
    float maxx = std::max(v1.x, std::max(v2.x, v3.x));
    float maxy = std::max(v1.y, std::max(v2.y, v3.y));
    float maxz = std::max(v1.z, std::max(v2.z, v3.z));
    return glm::vec3(maxx, maxy, maxz);
}

// loading a scene from obj file
// ---------------------------
class Scene
{
public:
    std::vector<Texture> textures_loaded;
    std::vector<Triangle> meshTriangles;
    std::vector<Triangle> lightTriangles;
    std::vector<Material> materials;
    std::vector<Vertex> vertices;
    std::vector<BVHNode> meshNodes;
    std::vector<BVHNode> lightNodes;
    std::vector<unsigned int> indices;
    glm::mat4 model;
    std::string directory;
    bool gammaCorrection;
    unsigned int VAO, VBO, EBO;
    unsigned int verticesTextureBuffer, meshTrianglesTextureBuffer, lightTrianglesTextureBuffer, materialsTextureBuffer, meshBVHTextureBuffer;
    unsigned int texture_array;
    const int MAXWIDTH = 4096, MAXHEIGHT = 4096, MAXNUM = 32;

    // constructor, expects a filepath to a 3D model.
    Scene(std::string const &path, unsigned int pFlags, glm::mat4 model = glm::mat4(1.0f), bool gamma = false) : model(model), gammaCorrection(gamma)
    {
        glGenTextures(1, &texture_array);
        glActiveTexture(GL_TEXTURE0 + texture_array);
        glBindTexture(GL_TEXTURE_2D_ARRAY, texture_array);
        glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_RGBA, MAXWIDTH, MAXHEIGHT, MAXNUM, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);

        loadScene(path, pFlags);

        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glGenerateMipmap(GL_TEXTURE_2D_ARRAY);
        glActiveTexture(0);

        buildBVHwithSAH(meshTriangles, meshNodes, 0, meshTriangles.size() - 1, 1);

        loadVertices();
        loadMeshTriangles();
        loadLightTriangles();
        loadMaterials();
        loadMeshBVHNodes();
    }

    ~Scene()
    {
        glDeleteTextures(1, &texture_array);
        glDeleteTextures(1, &meshBVHTextureBuffer);
        glDeleteTextures(1, &meshTrianglesTextureBuffer);
        glDeleteTextures(1, &verticesTextureBuffer);
        glDeleteVertexArrays(1, &VAO);
        glDeleteBuffers(1, &VBO);
        glDeleteBuffers(1, &EBO);
    }

    void loadShader(Shader &shader)
    {
        shader.use();
        shader.setInt("texturesArray", texture_array);
        shader.setInt("vertices", verticesTextureBuffer);
        shader.setInt("meshTriangles", meshTrianglesTextureBuffer);
        shader.setInt("lightTriangles", lightTrianglesTextureBuffer);
        shader.setInt("nlightTriangles", lightTriangles.size());
        shader.setInt("materials", materialsTextureBuffer);
        shader.setInt("meshBVHNodes", meshBVHTextureBuffer);
    }

    // draws the scene
    void draw()
    {
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);
    }
    
private:
    void loadScene(std::string const &path, unsigned int pFlags)
    {
        // read file via ASSIMP
        Assimp::Importer importer;
        const aiScene* scene = importer.ReadFile(path, pFlags);
        // check for errors
        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) { // if is Not Zero 
            std::cout << "ERROR::ASSIMP:: " << importer.GetErrorString() << std::endl;
            return;
        }
        // retrieve the directory path of the filepath
        directory = path.substr(0, path.find_last_of('/'));

        // process ASSIMP's root node recursively
        processNode(scene->mRootNode, scene);

        setupScene();
    }

    // load triangles into texture buffer
    void loadVertices()
    {
        GLuint tbo;
        glGenBuffers(1, &tbo);
        glBindBuffer(GL_TEXTURE_BUFFER, tbo);
        glBufferData(GL_TEXTURE_BUFFER, vertices.size() * sizeof(Vertex), &vertices[0], GL_STATIC_DRAW);
        glGenTextures(1, &verticesTextureBuffer);
        glActiveTexture(GL_TEXTURE0 + verticesTextureBuffer);
        glBindTexture(GL_TEXTURE_BUFFER, verticesTextureBuffer);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, tbo);
        glActiveTexture(0);
    }

    void loadMeshTriangles()
    {
        GLuint tbo;
        glGenBuffers(1, &tbo);
        glBindBuffer(GL_TEXTURE_BUFFER, tbo);
        glBufferData(GL_TEXTURE_BUFFER, meshTriangles.size() * sizeof(Triangle), &meshTriangles[0], GL_STATIC_DRAW);
        glGenTextures(1, &meshTrianglesTextureBuffer);
        glActiveTexture(GL_TEXTURE0 + meshTrianglesTextureBuffer);
        glBindTexture(GL_TEXTURE_BUFFER, meshTrianglesTextureBuffer);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, tbo);
        glActiveTexture(0);
    }

    void loadLightTriangles()
    {
        GLuint tbo;
        glGenBuffers(1, &tbo);
        glBindBuffer(GL_TEXTURE_BUFFER, tbo);
        glBufferData(GL_TEXTURE_BUFFER, lightTriangles.size() * sizeof(Triangle), &lightTriangles[0], GL_STATIC_DRAW);
        glGenTextures(1, &lightTrianglesTextureBuffer);
        glActiveTexture(GL_TEXTURE0 + lightTrianglesTextureBuffer);
        glBindTexture(GL_TEXTURE_BUFFER, lightTrianglesTextureBuffer);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, tbo);
        glActiveTexture(0);
    }

    void loadMaterials()
    {
        GLuint tbo;
        glGenBuffers(1, &tbo);
        glBindBuffer(GL_TEXTURE_BUFFER, tbo);
        glBufferData(GL_TEXTURE_BUFFER, materials.size() * sizeof(Material), &materials[0], GL_STATIC_DRAW);
        glGenTextures(1, &materialsTextureBuffer);
        glActiveTexture(GL_TEXTURE0 + materialsTextureBuffer);
        glBindTexture(GL_TEXTURE_BUFFER, materialsTextureBuffer);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, tbo);
        glActiveTexture(0);
    }

    void loadMeshBVHNodes()
    {
        GLuint tbo;
        glGenBuffers(1, &tbo);
        glBindBuffer(GL_TEXTURE_BUFFER, tbo);
        glBufferData(GL_TEXTURE_BUFFER, meshNodes.size() * sizeof(BVHNode), &meshNodes[0], GL_STATIC_DRAW);
        glGenTextures(1, &meshBVHTextureBuffer);
        glActiveTexture(GL_TEXTURE0 + meshBVHTextureBuffer);
        glBindTexture(GL_TEXTURE_BUFFER, meshBVHTextureBuffer);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, tbo);
        glActiveTexture(0);
    }

    // processes a node in a recursive fashion. Processes each individual mesh located at the node and repeats this process on its children nodes (if any).
    void processNode(aiNode *node, const aiScene *scene)
    {
        // process each mesh located at the current node
        for (unsigned int i = 0; i < node->mNumMeshes; i++) {
            // the node object only contains indices to index the actual objects in the scene. 
            // the scene contains all the data, node is just to keep stuff organized (like relations between nodes).
            aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
            processMesh(mesh, scene);
        }
        // after we've processed all of the meshes (if any) we then recursively process each of the children nodes
        for (unsigned int i = 0; i < node->mNumChildren; i++) {
            processNode(node->mChildren[i], scene);
        }
    }

    void processMesh(aiMesh *mesh, const aiScene *scene)
    {
        // process material
        aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
        Material mat;

        mat.ambient = loadMaterialTextures(material, aiTextureType_AMBIENT, AI_MATKEY_COLOR_AMBIENT);
        mat.diffuse = loadMaterialTextures(material, aiTextureType_DIFFUSE, AI_MATKEY_COLOR_DIFFUSE);
        mat.specular = loadMaterialTextures(material, aiTextureType_SPECULAR, AI_MATKEY_COLOR_SPECULAR);
        mat.emissive = loadMaterialTextures(material, aiTextureType_EMISSIVE, AI_MATKEY_COLOR_EMISSIVE);
        mat.shininess = loadMaterialTextures(material, aiTextureType_SHININESS, AI_MATKEY_SHININESS);
        mat.metallic = loadMaterialTextures(material, aiTextureType_METALNESS, AI_MATKEY_METALLIC_FACTOR);
        mat.refracti = loadMaterialTextures(material, aiTextureType_REFLECTION, AI_MATKEY_REFRACTI);
        mat.opacity = loadMaterialTextures(material, aiTextureType_OPACITY, AI_MATKEY_OPACITY);
        mat.transmission = loadMaterialTextures(material, aiTextureType_TRANSMISSION, AI_MATKEY_COLOR_TRANSPARENT);
        mat.anisotropy = loadMaterialTextures(material, aiTextureType_UNKNOWN, AI_MATKEY_ANISOTROPY_FACTOR);

        int matID = this->materials.size();
        this->materials.push_back(mat);

        glm::vec3 norm = loadMaterialTextures(material, aiTextureType_NORMALS);
        glm::vec3 bump = loadMaterialTextures(material, aiTextureType_HEIGHT);

        // walk through each of the mesh's faces (a face is a mesh its triangle) and retrieve the corresponding vertex indices.
        for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
            aiFace face = mesh->mFaces[i];
            Triangle t;
            // retrieve all indices of the face and store them in the indices vector
            for (unsigned int j = 0; j < face.mNumIndices; j ++) {
                t.v[j] = face.mIndices[j] + this->vertices.size();
                this->indices.push_back(face.mIndices[j] + this->vertices.size());
            }
            this->meshTriangles.push_back(t);
            if (mat.emissive.x > 0 || mat.emissive.y > 0 || mat.emissive.z > 0) {
                this->lightTriangles.push_back(t);
            }
        }

        // walk through each of the mesh's vertices
        for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
            Vertex vertex;

            glm::vec3 vector;
            glm::vec2 vec;
            // position
            vector.x = mesh->mVertices[i].x;
            vector.y = mesh->mVertices[i].y;
            vector.z = mesh->mVertices[i].z;
            vertex.pos = glm::vec3(this->model * glm::vec4(vector, 1));
            // normal
            if (norm.y < 0 && norm.z < 0) {
                vector.x = mesh->mNormals[i].x;
                vector.y = mesh->mNormals[i].y;
                vector.z = mesh->mNormals[i].z;
                vertex.normal = glm::mat3(this->model) * vector;
            }
            else {
                vertex.normal = norm;
            }
            // texcoord
            vec.x = mesh->mTextureCoords[0][i].x;
            vec.y = mesh->mTextureCoords[0][i].y;
            vertex.texCoords = vec;
            // matID
            vertex.materialID = matID;
            // tangent
            vector.x = mesh->mTangents[i].x;
            vector.y = mesh->mTangents[i].y;
            vector.z = mesh->mTangents[i].z;
            vertex.tangent = glm::mat3(this->model) * vector;
            // bitangent
            vector.x = mesh->mBitangents[i].x;
            vector.y = mesh->mBitangents[i].y;
            vector.z = mesh->mBitangents[i].z;
            vertex.bitangent = glm::mat3(this->model) * vector;

            this->vertices.push_back(vertex);
        }
    }

    // load material from file or rgb into texture
    glm::vec3 loadMaterialTextures(aiMaterial *material, aiTextureType type, const char *pKey = nullptr, unsigned int key_type = 0, unsigned int idx = 0)
    {
        if (material->GetTextureCount(type)) {
            aiString str;
            material->GetTexture(type, 0, &str); // just load the first texture
            // check if texture was loaded before and if so, continue to next iteration: skip loading a new texture
            for (unsigned int i = 0; i < textures_loaded.size(); i ++) {
                if (std::strcmp(textures_loaded[i].name.data(), str.C_Str()) == 0) {
                    return textures_loaded[i].id;
                }
            }
            Texture t;
            t.name = str.C_Str();
            t.id = TextureFromFile(str.C_Str());
            textures_loaded.push_back(t);
            return t.id;
        }
        else if (pKey != nullptr) {
            aiColor3D color;
            material->Get(pKey, key_type, idx, color);
            // std::cout << pKey << color.r << color.g << color.b << std::endl;
            return glm::vec3(color.r, color.g, color.b);
        }
        else {
            return glm::vec3(-1);
        }
    }

    // load texture from file
    glm::vec3 TextureFromFile(const char *path)
    {
        std::string filename = std::string(path);
        filename = directory + '/' + filename;

        int width, height, nrComponents;
        unsigned char *data = stbi_load(filename.c_str(), &width, &height, &nrComponents, 0);
        if (data  && width <= MAXWIDTH && height <= MAXHEIGHT) {
            GLenum format;
            if (nrComponents == 1)
                format = GL_RED;
            else if (nrComponents == 3)
                format = GL_RGB;
            else if (nrComponents == 4)
                format = GL_RGBA;

            glTexSubImage3D(GL_TEXTURE_2D_ARRAY, 0, 0, 0, textures_loaded.size(), width, height, 1, format, GL_UNSIGNED_BYTE, data);
            stbi_image_free(data);
            return glm::vec3(-2.0 - textures_loaded.size(), (float)width / MAXWIDTH, (float)height / MAXHEIGHT);
        }
        else {
            std::cout << "Texture failed to load at path: " << path << std::endl;
            stbi_image_free(data);
            exit(1);
        }
    }

    void setupScene()
    {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);

        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), &vertices[0], GL_STATIC_DRAW);  

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);

        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(6 * sizeof(float)));
        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(8 * sizeof(float)));

        glBindVertexArray(0);
    }

    // BVH with SAH
    int buildBVHwithSAH(std::vector<Triangle>& triangles, std::vector<BVHNode>& bvhNodes, int l, int r, int n) {
        if (l > r) return 0;

        BVHNode root;
        bvhNodes.push_back(root);
        int id = bvhNodes.size() - 1;
        bvhNodes[id].left = bvhNodes[id].right = bvhNodes[id].n = bvhNodes[id].index1 = bvhNodes[id].index2 = bvhNodes[id].parent = 0;
        bvhNodes[id].AA = this->vertices[triangles[l].v[0]].pos;
        bvhNodes[id].BB = this->vertices[triangles[l].v[0]].pos;

        // count AABB
        for (int i = l; i <= r; i++) {
            bvhNodes[id].AA = minVec3(bvhNodes[id].AA, this->vertices[triangles[i].v[0]].pos, this->vertices[triangles[i].v[1]].pos);
            bvhNodes[id].AA = minVec3(bvhNodes[id].AA, this->vertices[triangles[i].v[1]].pos, this->vertices[triangles[i].v[2]].pos);
            bvhNodes[id].BB = maxVec3(bvhNodes[id].BB, this->vertices[triangles[i].v[0]].pos, this->vertices[triangles[i].v[1]].pos);
            bvhNodes[id].BB = maxVec3(bvhNodes[id].BB, this->vertices[triangles[i].v[1]].pos, this->vertices[triangles[i].v[2]].pos);
        }

        // nTriangles <= n
        if ((r - l + 1) <= n) {
            bvhNodes[id].n = r - l + 1;
            bvhNodes[id].index1 = l;
            bvhNodes[id].index2 = r;
            return id;
        }

        // build
        int minCost = 0;
        int Axis = 0;
        int Split = (l + r) / 2;
        for (int axis = 0; axis < 3; axis ++) {
            // sort in x，y，z axis
            if (axis == 0) std::sort(&triangles[0] + l, &triangles[0] + r + 1, [this](const Triangle& t1, const Triangle& t2){return cmpx(this->vertices, t1, t2);});
            if (axis == 1) std::sort(&triangles[0] + l, &triangles[0] + r + 1, [this](const Triangle& t1, const Triangle& t2){return cmpy(this->vertices, t1, t2);});
            if (axis == 2) std::sort(&triangles[0] + l, &triangles[0] + r + 1, [this](const Triangle& t1, const Triangle& t2){return cmpz(this->vertices, t1, t2);});

            std::vector<glm::vec3> leftAA = {minVec3(this->vertices[triangles[l].v[0]].pos, this->vertices[triangles[l].v[1]].pos, this->vertices[triangles[l].v[2]].pos)};
            std::vector<glm::vec3> leftBB = {maxVec3(this->vertices[triangles[l].v[0]].pos, this->vertices[triangles[l].v[1]].pos, this->vertices[triangles[l].v[2]].pos)};
            for (int i = l + 1; i < r; i ++) {
                glm::vec3 AA = minVec3(leftAA[i - l - 1], this->vertices[triangles[i].v[0]].pos, this->vertices[triangles[i].v[1]].pos);
                AA = minVec3(AA, this->vertices[triangles[i].v[1]].pos, this->vertices[triangles[i].v[2]].pos);
                glm::vec3 BB = maxVec3(leftBB[i - l - 1], this->vertices[triangles[i].v[0]].pos, this->vertices[triangles[i].v[1]].pos);
                BB = maxVec3(BB, this->vertices[triangles[i].v[1]].pos, this->vertices[triangles[i].v[2]].pos);

                leftAA.push_back(AA);
                leftBB.push_back(BB);
            }

            std::vector<glm::vec3> rightAA = {minVec3(this->vertices[triangles[r].v[0]].pos, this->vertices[triangles[r].v[1]].pos, this->vertices[triangles[r].v[2]].pos)};
            std::vector<glm::vec3> rightBB = {maxVec3(this->vertices[triangles[r].v[0]].pos, this->vertices[triangles[r].v[1]].pos, this->vertices[triangles[r].v[2]].pos)};
            for (int i = r - 1; i > l; i --) {
                glm::vec3 AA = minVec3(rightAA[r - i - 1], this->vertices[triangles[i].v[0]].pos, this->vertices[triangles[i].v[1]].pos);
                AA = minVec3(AA, this->vertices[triangles[i].v[1]].pos, this->vertices[triangles[i].v[2]].pos);
                glm::vec3 BB = maxVec3(rightBB[r - i - 1], this->vertices[triangles[i].v[0]].pos, this->vertices[triangles[i].v[1]].pos);
                BB = maxVec3(BB, this->vertices[triangles[i].v[1]].pos, this->vertices[triangles[i].v[2]].pos);

                rightAA.push_back(AA);
                rightBB.push_back(BB);
            }

            // search for split
            float cost = 0;
            int split = l;
            for (int i = l; i < r; i ++) {
                float lenx = leftBB[i - l].x - leftAA[i - l].x;
                float leny = leftBB[i - l].y - leftAA[i - l].y;
                float lenz = leftBB[i - l].z - leftAA[i - l].z;
                float leftS = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny * lenz));
                float leftCost = leftS * (i - l + 1);

                lenx = rightBB[r - 1 - i].x - rightAA[r - 1 - i].x;
                leny = rightBB[r - 1 - i].y - rightAA[r - 1 - i].y;
                lenz = rightBB[r - 1 - i].z - rightAA[r - 1 - i].z;
                float rightS = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny * lenz));
                float rightCost = rightS * (r - i);

                float totalCost = leftCost + rightCost;
                if (cost == 0 || totalCost < cost) {
                    cost = totalCost;
                    split = i;
                }
            }
            if (minCost == 0 || cost < minCost) {
                minCost = cost;
                Axis = axis;
                Split = split;
            }
        }

        if (Axis == 0) std::sort(&triangles[0] + l, &triangles[0] + r + 1, [this](const Triangle& t1, const Triangle& t2){return cmpx(this->vertices, t1, t2);});
        if (Axis == 1) std::sort(&triangles[0] + l, &triangles[0] + r + 1, [this](const Triangle& t1, const Triangle& t2){return cmpy(this->vertices, t1, t2);});
        if (Axis == 2) std::sort(&triangles[0] + l, &triangles[0] + r + 1, [this](const Triangle& t1, const Triangle& t2){return cmpz(this->vertices, t1, t2);});

        int left  = buildBVHwithSAH(triangles, bvhNodes, l, Split, n);
        int right = buildBVHwithSAH(triangles, bvhNodes, Split + 1, r, n);

        bvhNodes[left].parent = id;
        bvhNodes[right].parent = id;
        bvhNodes[id].left = left;
        bvhNodes[id].right = right;

        return id;
    }
};
#endif