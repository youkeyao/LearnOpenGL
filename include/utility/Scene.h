#ifndef MODEL_H
#define MODEL_H

#include <glad/glad.h> 

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include <utility/stb_image.h>
#include <utility/Shader.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
using namespace std;

struct Texture {
    glm::vec3 id;
    string name;
};

struct Material
{
    glm::vec3 ambient;          // (id, width, height)
    glm::vec3 diffuse;          // (id, width, height)
    glm::vec3 specular;         // (id, width, height)
    glm::vec3 emissive;         // (id, width, height)
    glm::vec3 shininess;        // (id, width, height)
    glm::vec3 metallic;         // (id, width, height)
    glm::vec3 refracti;         // (id, width, height)
    glm::vec3 opacity;          // (id, width, height)
    glm::vec3 transmission;     // (id, width, height)
    glm::vec3 anisotropy;       // (id, width, height)
};

struct Triangle {
    glm::vec3 p[3];
    glm::vec3 n[3];
    glm::vec2 texCoords[3];
    glm::mat3 TBN[3];
    glm::vec3 normal;           // (id, width, height)
    glm::vec3 height;           // (id, width, height)
    Material material;
};

class Scene
{
public:
    vector<Texture> textures_loaded;
    vector<Triangle> triangles;
    vector<glm::vec3> vertices;
    vector<unsigned int> indices;
    glm::mat4 model;
    string directory;
    bool gammaCorrection;
    unsigned int VAO, VBO, EBO, trianglesTextureBuffer, texture_array;
    const int MAXWIDTH = 4096, MAXHEIGHT = 4096, MAXNUM = 32;

    // constructor, expects a filepath to a 3D model.
    Scene(string const &path, unsigned int pFlags, glm::mat4 model = glm::mat4(1.0f), bool gamma = false) : model(model), gammaCorrection(gamma)
    {
        glGenTextures(1, &texture_array);
        glBindTexture(GL_TEXTURE_2D_ARRAY, texture_array);
        glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_RGBA, MAXWIDTH, MAXHEIGHT, MAXNUM, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);

        loadScene(path, pFlags);

        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glGenerateMipmap(GL_TEXTURE_2D_ARRAY);

        loadTriangles();
    }

    // draws the scene
    void Draw(Shader &shader)
    {
        // active proper texture unit before binding
        
        // now set the sampler to the correct texture unit
        glActiveTexture(GL_TEXTURE0 + trianglesTextureBuffer - 1);
        glBindTexture(GL_TEXTURE_BUFFER, trianglesTextureBuffer);
        shader.setInt("triangles", trianglesTextureBuffer - 1);
        shader.setInt("nTriangles", triangles.size());

        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        glActiveTexture(0);
    }
    
private:
    void loadScene(string const &path, unsigned int pFlags)
    {
        // read file via ASSIMP
        Assimp::Importer importer;
        const aiScene* scene = importer.ReadFile(path, pFlags);
        // check for errors
        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) { // if is Not Zero 
            cout << "ERROR::ASSIMP:: " << importer.GetErrorString() << endl;
            return;
        }
        // retrieve the directory path of the filepath
        directory = path.substr(0, path.find_last_of('/'));

        // process ASSIMP's root node recursively
        processNode(scene->mRootNode, scene);

        setupScene();
    }

    // load triangles into shader
    void loadTriangles()
    {
        GLuint tbo0;
        glGenBuffers(1, &tbo0);
        glBindBuffer(GL_TEXTURE_BUFFER, tbo0);
        glBufferData(GL_TEXTURE_BUFFER, triangles.size() * sizeof(Triangle), &triangles[0], GL_STATIC_DRAW);
        glGenTextures(1, &trianglesTextureBuffer);
        glActiveTexture(GL_TEXTURE0 + trianglesTextureBuffer - 1);
        glBindTexture(GL_TEXTURE_BUFFER, trianglesTextureBuffer);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, tbo0);
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
        mat.transmission = loadMaterialTextures(material, aiTextureType_TRANSMISSION, AI_MATKEY_TRANSMISSION_FACTOR);
        mat.anisotropy = loadMaterialTextures(material, aiTextureType_UNKNOWN, AI_MATKEY_ANISOTROPY_FACTOR);

        glm::vec3 norm = loadMaterialTextures(material, aiTextureType_NORMALS);
        glm::vec3 bump = loadMaterialTextures(material, aiTextureType_HEIGHT);

        // walk through each of the mesh's faces (a face is a mesh its triangle) and retrieve the corresponding vertex indices.
        for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
            aiFace face = mesh->mFaces[i];
            Triangle t;
            t.material = mat;
            t.normal = norm;
            t.height = bump;
            // retrieve all indices of the face and store them in the indices vector
            for (unsigned int j = 0; j < face.mNumIndices; j ++) {
                glm::vec3 vector, tangent, bitangent;
                glm::vec2 vec;
                // position
                vector.x = mesh->mVertices[face.mIndices[j]].x;
                vector.y = mesh->mVertices[face.mIndices[j]].y;
                vector.z = mesh->mVertices[face.mIndices[j]].z;
                t.p[j] = glm::vec3(this->model * glm::vec4(vector, 0.0));
                // normal
                vector.x = mesh->mNormals[face.mIndices[j]].x;
                vector.y = mesh->mNormals[face.mIndices[j]].y;
                vector.z = mesh->mNormals[face.mIndices[j]].z;
                t.n[j] = glm::normalize(glm::transpose(glm::inverse(glm::mat3(this->model))) * vector);
                // texture coordinates
                vec.x = mesh->mTextureCoords[0][face.mIndices[j]].x;
                vec.y = mesh->mTextureCoords[0][face.mIndices[j]].y;
                t.texCoords[j] = vec;
                // tangent
                tangent.x = mesh->mTangents[face.mIndices[j]].x;
                tangent.y = mesh->mTangents[face.mIndices[j]].y;
                tangent.z = mesh->mTangents[face.mIndices[j]].z;
                // bitangent
                bitangent.x = mesh->mBitangents[face.mIndices[j]].x;
                bitangent.y = mesh->mBitangents[face.mIndices[j]].y;
                bitangent.z = mesh->mBitangents[face.mIndices[j]].z;

                t.TBN[j] = calcTBN(tangent, bitangent, t.n[j]);

                this->indices.push_back(face.mIndices[j] + this->vertices.size());
            }
            triangles.push_back(t);
        }

        // walk through each of the mesh's vertices
        for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
            glm::vec3 vector;
            vector.x = mesh->mVertices[i].x;
            vector.y = mesh->mVertices[i].y;
            vector.z = mesh->mVertices[i].z;

            vertices.push_back(glm::vec3(this->model * glm::vec4(vector, 0.0)));
        }
    }

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
            Texture t;
            t.name = color2string(color);
            t.id = TextureFromRGB(color.r, color.g, color.b);
            textures_loaded.push_back(t);
            return t.id;
        }
        else {
            return glm::vec3(-1, 0, 0);
        }
    }

    void setupScene()
    {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);

        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), &vertices[0], GL_STATIC_DRAW);  

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);

        glEnableVertexAttribArray(0);	
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);

        glBindVertexArray(0);
    }

    glm::mat3 calcTBN(glm::vec3 tangent, glm::vec3 bitangent, glm::vec3 normal)
    {
        glm::vec3 T = glm::normalize(glm::transpose(glm::inverse(glm::mat3(this->model))) * tangent);
        glm::vec3 B = glm::normalize(glm::transpose(glm::inverse(glm::mat3(this->model))) * bitangent);
        glm::vec3 N = glm::normalize(normal);
        return glm::mat3(T, B, N);
    }

    glm::vec3 TextureFromRGB(float r, float g, float b)
    {
        unsigned char data[] = { (unsigned char)(r * 255), (unsigned char)(g * 255), (unsigned char)(b * 255), 255 };
        glTexSubImage3D(GL_TEXTURE_2D_ARRAY, 0, 0, 0, textures_loaded.size(), 1, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, data);
        return glm::vec3(textures_loaded.size(), 1 / (float)MAXWIDTH, 1 / (float)MAXHEIGHT);
    }

    glm::vec3 TextureFromFile(const char *path)
    {
        string filename = string(path);
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
            return glm::vec3(textures_loaded.size(), (float)width / MAXWIDTH, (float)height / MAXHEIGHT);
        }
        else {
            std::cout << "Texture failed to load at path: " << path << std::endl;
            stbi_image_free(data);
            exit(1);
        }
    }

    string color2string(aiColor3D color)
    {
        string s;
        char buf[30];
        sprintf(buf, "%.4f %.4f %.4f", color.r, color.g, color.b);
        s.assign(buf);
        return s;
    }
};
#endif