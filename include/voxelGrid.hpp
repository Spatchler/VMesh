#pragma once

#include "mesh.hpp"

#include <future>
#include <sstream>
#include <fstream>
#include <iostream>
#include <functional>
#include <bitset>

namespace VMesh {
  constexpr float& v3index(glm::tvec3<float>& v, uint8_t i) {
    if (i == 0) return v.x;
    if (i == 1) return v.y;
    if (i == 2) return v.z;
    return v.x;
  }

  typedef std::function<void(float,float,float)> insertFunc_t;

  class VoxelGrid {
  public:
    VoxelGrid(uint pResolution, glm::vec3 pOrigin = glm::vec3(0,0,0));

    void setOrigin(glm::vec3 pOrigin = glm::vec3(0,0,0));
    glm::vec3 getOrigin();

    void voxelizeMesh(Mesh& pMesh, uint* pTrisComplete);
    void DDAvoxelizeMesh(Mesh& pMesh, uint* pTrisComplete , const insertFunc_t& pInsertFunc = [](...){});

    void writeToFile(const std::string& pPath);
    void writeToFileCompressed(const std::string& pPath, uint64_t* pVoxelsComplete);

    void loadFromFile(const std::string& pPath);
    void loadFromFileCompressed(const std::string& pPath);

    void setLogStream(std::ostream* pStream, std::mutex* pMutex = NULL);

    int insert(const glm::uvec3& pPos, const insertFunc_t& pInsertFunc = [](...){});

    void voxelizeTriangle(std::array<glm::vec3, 3> pPoints);

    uint64_t getVoxelCount();
    uint64_t getVolume();
    uint getResolution();
    float getMaxDepth();
    bool queryVoxel(const glm::uvec3& pPos);
    // const std::vector<std::vector<std::vector<bool>>>& getVoxelData();
    const std::vector<char>& getVoxelDataBits();
    std::vector<uint64_t> generateCompressedVoxelData(uint64_t* pVoxelsComplete);

    std::mutex mDefaultLogMutex;
    std::stringstream mDefaultLogStream;
  protected:
    void init();

    bool triangleBoxIntersect(const std::array<glm::vec3, 3>& pTriPoints, const glm::vec3& pBoxOrigin);
    bool planeBoxOverlap(glm::vec3& pNormal, glm::vec3& pVert);
    void findMinMax(float& a, float& b, float& c, float& pMin, float& pMax);
    bool axistestX01(const glm::vec3& pTranslatedTriPoint0, const glm::vec3& pTranslatedTriPoint2, float a, float b, float fa, float fb, float& pMin, float& pMax);
    bool axistestX2(const glm::vec3& pTranslatedTriPoint0, const glm::vec3& pTranslatedTriPoint1, float a, float b, float fa, float fb, float& pMin, float& pMax);
    bool axistestY02(const glm::vec3& pTranslatedTriPoint0, const glm::vec3& pTranslatedTriPoint2, float a, float b, float fa, float fb, float& pMin, float& pMax);
    bool axistestY1(const glm::vec3& pTranslatedTriPoint0, const glm::vec3& pTranslatedTriPoint1, float a, float b, float fa, float fb, float& pMin, float& pMax);
    bool axistestZ12(const glm::vec3& pTranslatedTriPoint1, const glm::vec3& pTranslatedTriPoint2, float a, float b, float fa, float fb, float& pMin, float& pMax);
    bool axistestZ0(const glm::vec3& pTranslatedTriPoint1, const glm::vec3& pTranslatedTriPoint2, float a, float b, float fa, float fb, float& pMin, float& pMax);

    void openFileWrite(std::ofstream& pFout, const std::string& pPath);
    void writeMetaData(std::ofstream& pFout);

    static glm::vec3 toVec3(float a, float b, float c, uint8_t pDominantAxisIndex);

    void drawLine2(uint8_t pDominantAxisIndex, float pDominantAxisValue, const glm::vec2& pStart, const glm::vec2& pEnd, const glm::vec2& pDir, const glm::vec2& pDirInv);

    std::ostream* mLogStream;
    std::mutex* mLogMutex;

    // std::vector<std::vector<std::vector<bool>>> mVoxelGrid;
    std::vector<char> mVoxelData;
    
    glm::vec3 mOrigin;
    uint mResolution = 0;
    uint64_t mVolume = 0, mVoxelCount = 0;
    float mMaxDepth = 0;
  };
}

