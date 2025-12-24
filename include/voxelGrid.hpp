#pragma once

#include "mesh.hpp"

#include <future>
#include <sstream>
#include <fstream>
#include <iostream>
#include <functional>

namespace VMesh {
  float v3index(glm::tvec3<float> v, uint8_t i);

  typedef std::function<void(float,float,float)> insertFunc_t;

  class VoxelGrid {
  public:
    VoxelGrid(uint pResolution);

    void voxelizeMesh(Mesh& pMesh, uint* pTrisComplete = NULL, const insertFunc_t& pInsertFunc = [](...){});

    void writeToFile(const std::string& pPath);
    void writeToFileCompressed(const std::string& pPath);

    void setLogStream(std::ostream* pStream, std::mutex* pMutex = NULL);

    int insert(const glm::vec3& pPos, const insertFunc_t& pInsertFunc = [](...){});

    uint getVoxelCount();
    const std::vector<std::vector<std::vector<bool>>>& getVoxelData();
    std::vector<uint> generateCompressedVoxelData();

    std::mutex mDefaultLogMutex;
  protected:
    void openFileWrite(std::ofstream& pFout, const std::string& pPath);
    void writeMetaData(std::ofstream& pFout);

    static glm::vec3 toVec3(float a, float b, float c, uint8_t pDominantAxisIndex);

    void drawLine2(uint8_t pDominantAxisIndex, float pDominantAxisValue, const glm::vec2& pStart, const glm::vec2& pEnd, const glm::vec2& pDir, const glm::vec2& pDirInv, const insertFunc_t& pInsertFunc);

    std::ostream* mLogStream;
    std::mutex* mLogMutex;

    std::vector<std::vector<std::vector<bool>>> mVoxelGrid;
    
    uint mResolution, mVoxelCount = 0;
  };
}

