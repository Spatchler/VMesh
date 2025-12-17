#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
// #include <glm/gtc/type_ptr.hpp>

#include "timer.hpp"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include <string>
#include <functional>
#include <sstream>
#include <iostream>
#include <future>

namespace VMesh {
  float v3index(glm::tvec3<float> v, uint8_t i);

  typedef std::function<void(float,float,float)> insertFunc_t;

  class Model {
  public:
    Model() = default;
    Model(const std::string& pPath, const uint& pResolution, const insertFunc_t& pInsertFunc = [](...){});

    void load(const std::string& pPath, const uint& pResolution, const insertFunc_t& pInsertFunc = [](...){});

    void loadMeshData(const std::string& pPath);

    void transformMeshVertices(const glm::mat4& pM);

    uint getTriCount();
    uint getVoxelCount();

    float getProgress();

    void setLogStream(std::ostream* pStream, std::mutex* pMutex = NULL);

    int insert(const glm::vec3& pPos, const insertFunc_t& pInsertFunc = [](...){});

    void generateVoxelData(uint pResolution, const insertFunc_t& pInsertFunc = [](...){});

    const std::vector<glm::vec3>& getMeshVertices();
    const std::vector<uint>& getMeshIndices();
    const std::vector<std::vector<std::vector<bool>>>& getVoxelData();
    std::vector<uint> generateCompressedVoxelData();

    std::mutex mDefaultLogMutex;
  protected:
    void processNode(aiNode* pNode, const aiScene* pScene);
    void processMesh(aiMesh* pMesh, const aiScene* pScene);

    static glm::vec3 toVec3(float a, float b, float c, uint8_t pDominantAxisIndex);

    void drawLine2(uint8_t pDominantAxisIndex, float pDominantAxisValue, const glm::vec2& pStart, const glm::vec2& pEnd, const glm::vec2& pDir, const glm::vec2& pDirInv, const insertFunc_t& pInsertFunc);

    std::ostream* mLogStream;
    std::mutex* mLogMutex;

    std::vector<glm::vec3> mMeshVertices;
    std::vector<uint> mMeshIndices;

    std::vector<std::vector<std::vector<bool>>> mVoxelGrid;

    float mTriCountInv;
    uint mResolution, mVoxelCount, mTriCount, mTrisComplete;
  };
}

