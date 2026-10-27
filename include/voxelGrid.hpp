#pragma once

#include "model.hpp"
#include "texture.hpp"
#include "palette.hpp"

#include <glm/gtc/type_ptr.hpp>

#include <cassert>

#include <filesystem>
#include <future>
#include <sstream>
#include <fstream>
#include <iostream>
#include <functional>
#include <bitset>

namespace VMesh {
  struct Triangle {
    std::array<Vertex, 3> verts;
    glm::vec3 normal;
    glm::vec3 absNormal;
    glm::mat3 barycentricMat;
  };

  struct VoxelGridMetadata {
    virtual void loadFromFile(std::ifstream& fin);
    virtual void print() noexcept;
    std::array<uint8_t, 3> fileVersion;
    uint32_t resolution;
    uint8_t paletteSize;
    uint64_t voxelCount;
  };

  struct CompressedVoxelGridMetadata: public VoxelGridMetadata {
    void loadFromFile(std::ifstream& fin) override;
    void print() noexcept override;
    uint32_t countsLen;
  };

  class VoxelGrid {
  public:
    VoxelGrid(uint pResolution, glm::vec3 pOrigin = glm::vec3(0,0,0));

    void setOrigin(glm::vec3 pOrigin = glm::vec3(0,0,0));
    glm::vec3 getOrigin();

    void clear();

    void IntersectVoxelizeMesh(Mesh& pMesh, uint* pTrisComplete);
    void IntersectVoxelizeModel(Model& pMesh, uint* pTrisComplete);
    void DDAvoxelizeMesh(Mesh& pMesh, uint* pTrisComplete, float pAddColourDistance, bool pCreatePalette, Texture* pTex = NULL);
    void DDAvoxelizeModel(Model& pModel, uint* pTrisComplete, bool pColoured, bool pCreatePalette, float pAddColourDistance);

    void writeToFile(std::string pPath);
    void writeToFileCompressed(std::string pPath, uint64_t* pVoxelsComplete);

    void loadFromFile(const std::string& pPath);
    void loadFromFileCompressed(const std::string& pPath);
    void loadFromVoxFile(const std::filesystem::path& pPath);

    void setLogStream(std::ostream* pStream, std::mutex* pMutex = NULL);

    void setUseMatrixBarycentricCompute(bool pState);

    void insert(uint64_t pIndex, uint8_t pCol);
    void insert(const glm::uvec3& pPos, uint8_t pCol);

    void DDAvoxelizeTriangle(std::array<Vertex, 3> pVerts, float pAddColourDistance, bool pCreatePalette, Texture* pTex = NULL);

    bool isOctreeOctantAllSame(const uint64_t& pIndex, uint pSize);
    bool isOctree2x2AllSame(const uint64_t& pIndex, uint8_t pData) noexcept;

    uint64_t getVoxelCount();
    uint64_t getVolume();
    uint getResolution();
    float getMaxDepth();
    bool queryVoxelPresence(uint64_t pIndex);
    bool queryVoxelPresence(const glm::uvec3& pPos);
    uint8_t queryVoxelData(uint64_t pIndex);
    uint8_t queryVoxelData(const glm::uvec3& pPos);
    const std::vector<uint8_t>& getVoxelData();
    std::pair<std::vector<uint8_t>, std::vector<uint64_t>> generateCompressedVoxelData(uint64_t* pVoxelsComplete);

    static uint64_t zorder(const glm::uvec3& pPos);

    Palette mPalette;

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

    void writeMetaData(std::ofstream& pFout);
    void readMetaData(std::ifstream& pFin);

    // Point and triangle vertices all in world space
    glm::vec3 computeBarycentric(const glm::vec3& pPoint, const Triangle& pTri);
    glm::vec2 computeTexCoord(const glm::vec3& pPoint, const Triangle& pTri);
    uint8_t sample(const glm::vec3& pPoint, const Triangle& pTri, bool pCreatePalette, float pAddColourDistance, Texture* pTex);

    static glm::vec3 toVec3(float a, float b, float c, uint8_t pDominantAxisIndex);

    void drawLine2(uint8_t pDominantAxisIndex, float pDominantAxisValue, const glm::vec2& pStart, const glm::vec2& pEnd, const glm::vec2& pDir, const glm::vec2& pDirInv, bool pCreatePalette, float pAddColourDistance, Texture* pTex, const Triangle& pTri);

    std::ostream* mLogStream;
    std::mutex* mLogMutex;

    bool mIsUseMatrixBarycentricComp { false };

    std::vector<uint8_t> mVoxelData;
    
    glm::vec3 mOrigin;
    uint mResolution = 0, mResolution2 = 0;
    uint64_t mVolume = 0, mVoxelCount = 0;
    float mMaxDepth = 0;
  };
}

