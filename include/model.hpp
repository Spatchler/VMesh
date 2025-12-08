#pragma once

#include "helpers.hpp"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include <string>
#include <functional>

namespace VMesh {
  class Model {
  public:
    Model(const std::string& pPath, uint pResolution);

    std::function<void(glm::vec3)> insertCallback;
  protected:
    void loadModel(const std::string& pPath);
    void processNode(aiNode* pNode, const aiScene* pScene);
    void processMesh(aiMesh* pMesh, const aiScene* pScene);

    static glm::vec3 toVec3(float a, float b, float c, uint8_t pDominantAxisIndex);

    void drawLine2(uint8_t pDominantAxisIndex, float pDominantAxisValue, const glm::vec2& pStart, const glm::vec2& pEnd, const glm::vec2& pDir, const glm::vec2& pDirInv);

    std::vector<glm::vec3> mMeshVertices;
    std::vector<uint> mMeshIndices;

    std::vector<std::vector<std::vector<bool>>> mVoxelGrid;
    uint mResolution;

    uint mVoxelCount;
    uint mTriangleCount;
  };
}

