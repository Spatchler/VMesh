#pragma once

#include "mesh.hpp"

#include <future>
#include <sstream>
#include <fstream>
#include <iostream>
#include <functional>
#include <bitset>

namespace VMesh {
  constexpr float v3index(glm::tvec3<float> v, uint8_t i) {
    if (i == 0)
      return v.x;
    if (i == 1)
      return v.y;
    if (i == 2)
      return v.z;
    return 0.f;
  }

  typedef std::function<void(float,float,float)> insertFunc_t;

  class VoxelGrid {
  public:
    VoxelGrid(uint pResolution);

    void voxelizeMesh(Mesh& pMesh, uint* pTrisComplete , const insertFunc_t& pInsertFunc = [](...){});

    void writeToFile(const std::string& pPath);
    void writeToFileCompressed(const std::string& pPath, uint* pVoxelsComplete);

    void loadFromFile(const std::string& pPath);
    void loadFromFileCompressed(const std::string& pPath);

    void setLogStream(std::ostream* pStream, std::mutex* pMutex = NULL);

    int insert(const glm::tvec3<uint>& pPos, const insertFunc_t& pInsertFunc = [](...){});

    uint getVoxelCount();
    uint getVolume();
    uint getResolution();
    float getMaxDepth();
    const std::vector<std::vector<std::vector<bool>>>& getVoxelData();
    const std::vector<char>& getVoxelDataBits();
    std::vector<uint> generateCompressedVoxelData(uint* pVoxelsComplete);

    std::mutex mDefaultLogMutex;
  protected:
    void init();

    void openFileWrite(std::ofstream& pFout, const std::string& pPath);
    void writeMetaData(std::ofstream& pFout);

    static glm::vec3 toVec3(float a, float b, float c, uint8_t pDominantAxisIndex);

    void drawLine2(uint8_t pDominantAxisIndex, float pDominantAxisValue, const glm::vec2& pStart, const glm::vec2& pEnd, const glm::vec2& pDir, const glm::vec2& pDirInv, const insertFunc_t& pInsertFunc);

    std::ostream* mLogStream;
    std::mutex* mLogMutex;

    std::vector<std::vector<std::vector<bool>>> mVoxelGrid;
    std::vector<char> mVoxelData;
    
    uint mResolution, mVoxelCount, mVolume = 0;
    float mMaxDepth = 0;
  };
}

