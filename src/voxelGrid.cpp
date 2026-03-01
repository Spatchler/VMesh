#include "voxelGrid.hpp"

using namespace VMesh;

VoxelGrid::VoxelGrid(uint pResolution, glm::vec3 pOrigin)
:mResolution(pResolution), mOrigin(pOrigin) {
  mLogStream = &mDefaultLogStream;
  mLogMutex = &mDefaultLogMutex;
  init();
}

void VoxelGrid::setOrigin(glm::vec3 pOrigin) {
  mOrigin = pOrigin;
}

glm::vec3 VoxelGrid::getOrigin() {
  return mOrigin;
}

void VoxelGrid::voxelizeMesh(Mesh& pMesh, uint* pTrisComplete) {
  const std::vector<glm::vec3>& verts = pMesh.getVertices();
  const std::vector<uint>& indices = pMesh.getIndices();

  std::array<glm::vec3, 3> points;

  for (uint i = 0; i < pMesh.getTriCount(); ++i, ++*pTrisComplete) { // Every triangle
    // Add points to array
    for (uint j = 0; j < 3; ++j) points[j] = verts[indices[i*3+j]] - mOrigin;

    glm::vec3 minPoint = glm::min(points[0], glm::min(points[1], points[2]));
    glm::vec3 maxPoint = glm::max(points[0], glm::max(points[1], points[2]));

    minPoint = glm::floor(minPoint);
    maxPoint = glm::ceil(maxPoint);

    minPoint = glm::max(minPoint, glm::vec3(0));
    maxPoint = glm::min(maxPoint, glm::vec3(mResolution));

    for (glm::uvec3 pos = minPoint; pos.z <= maxPoint.z; ++pos.z) for (pos.y = minPoint.y; pos.y <= maxPoint.y; ++pos.y) for (pos.x = minPoint.x; pos.x <= maxPoint.x; ++pos.x) // Every block
      // if (!queryVoxel(pos)) if (triangleBoxIntersect(points, glm::vec3(pos) + 0.5f)) insert(pos);
      if (triangleBoxIntersect(points, glm::vec3(pos) + 0.5f)) insert(pos);
  }
}

bool VoxelGrid::triangleBoxIntersect(const std::array<glm::vec3, 3>& pTriPoints, const glm::vec3& pBoxOrigin) {
  // Translate triangle so box is center o
  glm::vec3 translatedTriPoint0 = pTriPoints[0] - pBoxOrigin;
  glm::vec3 translatedTriPoint1 = pTriPoints[1] - pBoxOrigin;
  glm::vec3 translatedTriPoint2 = pTriPoints[2] - pBoxOrigin;
  // Compute triangle edges
  glm::vec3 e0 = pTriPoints[1] - pTriPoints[0];
  glm::vec3 e1 = pTriPoints[2] - pTriPoints[1];
  glm::vec3 e2 = pTriPoints[0] - pTriPoints[2];
  glm::vec3 fe = glm::abs(e0);

  float min, max;

  if (axistestX01(translatedTriPoint0, translatedTriPoint2, e0.z, e0.y, fe.z, fe.y, min, max)) return false;
  if (axistestY02(translatedTriPoint0, translatedTriPoint2, e0.z, e0.x, fe.z, fe.x, min, max)) return false;
  if (axistestZ12(translatedTriPoint1, translatedTriPoint2, e0.y, e0.x, fe.y, fe.x, min, max)) return false;

  fe = glm::abs(e1);

  if (axistestX01(translatedTriPoint0, translatedTriPoint2, e1.z, e1.y, fe.z, fe.y, min, max)) return false;
  if (axistestY02(translatedTriPoint0, translatedTriPoint2, e1.z, e1.x, fe.z, fe.x, min, max)) return false;
  if (axistestZ0(translatedTriPoint0, translatedTriPoint1, e1.y, e1.x, fe.y, fe.x, min, max)) return false;

  fe = glm::abs(e2);

  if (axistestX2(translatedTriPoint0, translatedTriPoint1, e2.z, e2.y, fe.z, fe.y, min, max)) return false;
  if (axistestY1(translatedTriPoint0, translatedTriPoint1, e2.z, e2.x, fe.z, fe.x, min, max)) return false;
  if (axistestZ12(translatedTriPoint1, translatedTriPoint2, e2.y, e2.x, fe.y, fe.x, min, max)) return false;

  // Test x
  // findMinMax(translatedTriPoint0.x, translatedTriPoint1.x, translatedTriPoint2.x, min, max);
  // if (min > 0.5f || max < -0.5f) return false;
  // // Text y
  // findMinMax(translatedTriPoint0.y, translatedTriPoint1.y, translatedTriPoint2.y, min, max);
  // if (min > 0.5f || max < -0.5f) return false;
  // // Test z
  // findMinMax(translatedTriPoint0.z, translatedTriPoint1.z, translatedTriPoint2.z, min, max);
  // if (min > 0.5f || max < -0.5f) return false;

  // Test if the box intersects the plane of the triangle
  glm::vec3 triNormal = glm::cross(e0, e1);
  if (!planeBoxOverlap(triNormal, translatedTriPoint0)) return false;

  return true; // Box and triangle overlap
}

bool VoxelGrid::planeBoxOverlap(glm::vec3& pNormal, glm::vec3& pVert) {
  glm::vec3 vmin, vmax;
  float v;
  for (uint8_t i = 0; i <= 2; ++i) {
    v = v3index(pVert, i);
    if (v3index(pNormal, i) > 0.0f) {
      v3index(vmin, i) = -0.5f - v;
      v3index(vmax, i) =  0.5f - v;
    }
    else {
      v3index(vmin, i) =  0.5f - v;
      v3index(vmax, i) = -0.5f - v;
    }
  }
  if (glm::dot(pNormal, vmin) >  0.0f) return false;
  if (glm::dot(pNormal, vmax) >= 0.0f) return true;

  return false;
}

void VoxelGrid::findMinMax(float& a, float& b, float& c, float& pMin, float& pMax) {
  pMin = pMax = a;
  if (b < pMin) pMin = b;
  if (b > pMax) pMax = b;
  if (c < pMin) pMin = c;
  if (c > pMax) pMax = c;
}

bool VoxelGrid::axistestX01(const glm::vec3& pTranslatedTriPoint0, const glm::vec3& pTranslatedTriPoint2, float a, float b, float fa, float fb, float& pMin, float& pMax) {
  float p0 = a*pTranslatedTriPoint0.y - b*pTranslatedTriPoint0.z;
  float p2 = a*pTranslatedTriPoint2.y - b*pTranslatedTriPoint2.z;
  if (p0 < p2) { pMin=p0; pMax=p2; } else { pMin=p2; pMax=p0; }
  float rad = fa * 0.5f + fb * 0.5f;
  if (pMin > rad || pMax < -rad) return true;
  return false;
}

bool VoxelGrid::axistestX2(const glm::vec3& pTranslatedTriPoint0, const glm::vec3& pTranslatedTriPoint1, float a, float b, float fa, float fb, float& pMin, float& pMax) {
  float p0 = a*pTranslatedTriPoint0.y - b*pTranslatedTriPoint0.z;
  float p1 = a*pTranslatedTriPoint1.y - b*pTranslatedTriPoint1.z;
  if (p0 < p1) { pMin=p0; pMax=p1; } else { pMin=p1; pMax=p0; }
  float rad = fa * 0.5f + fb * 0.5f;
  if (pMin > rad || pMax < -rad) return true;
  return false;
}

bool VoxelGrid::axistestY02(const glm::vec3& pTranslatedTriPoint0, const glm::vec3& pTranslatedTriPoint2, float a, float b, float fa, float fb, float& pMin, float& pMax) {
  float p0 = -a*pTranslatedTriPoint0.x + b*pTranslatedTriPoint0.z;
  float p2 = -a*pTranslatedTriPoint2.x + b*pTranslatedTriPoint2.z;
  if (p0 < p2) { pMin=p0; pMax=p2; } else { pMin=p2; pMax=p0; }
  float rad = fa * 0.5f + fb * 0.5f;
  if (pMin > rad || pMax < -rad) return true;
  return false;
}

bool VoxelGrid::axistestY1(const glm::vec3& pTranslatedTriPoint0, const glm::vec3& pTranslatedTriPoint1, float a, float b, float fa, float fb, float& pMin, float& pMax) {
  float p0 = -a*pTranslatedTriPoint0.x + b*pTranslatedTriPoint0.z;
  float p1 = -a*pTranslatedTriPoint1.x + b*pTranslatedTriPoint1.z;
  if (p0 < p1) { pMin=p0; pMax=p1; } else { pMin=p1; pMax=p0; }
  float rad = fa * 0.5f + fb * 0.5f;
  if (pMin > rad || pMax < -rad) return true;
  return false;
}

bool VoxelGrid::axistestZ12(const glm::vec3& pTranslatedTriPoint1, const glm::vec3& pTranslatedTriPoint2, float a, float b, float fa, float fb, float& pMin, float& pMax) {
  float p1 = a*pTranslatedTriPoint1.x - b*pTranslatedTriPoint1.y;
  float p2 = a*pTranslatedTriPoint2.x - b*pTranslatedTriPoint2.y;
  if (p2 < p1) { pMin=p2; pMax=p1; } else { pMin=p1; pMax=p2; }
  float rad = fa * 0.5f + fb * 0.5f;
  if (pMin > rad || pMax < -rad) return true;
  return false;
}

bool VoxelGrid::axistestZ0(const glm::vec3& pTranslatedTriPoint0, const glm::vec3& pTranslatedTriPoint1, float a, float b, float fa, float fb, float& pMin, float& pMax) {
  float p0 = a*pTranslatedTriPoint0.x - b*pTranslatedTriPoint0.y;
  float p1 = a*pTranslatedTriPoint1.x - b*pTranslatedTriPoint1.y;
  if (p0 < p1) { pMin=p0; pMax=p1; } else { pMin=p1; pMax=p0; }
  float rad = fa * 0.5f + fb * 0.5f;
  if (pMin > rad || pMax < -rad) return true;
  return false;
}

void VoxelGrid::DDAvoxelizeMesh(Mesh& pMesh, uint* pTrisComplete, const insertFunc_t& pInsertFunc) {
  const std::vector<glm::vec3>& verts = pMesh.getVertices();
  const std::vector<uint>& indices = pMesh.getIndices();

  for (uint i = 0; i < pMesh.getTriCount(); ++i, ++*pTrisComplete) { // Every triangle
    std::array<glm::vec3, 3> points;

    bool xp = true, xn = true;
    bool yp = true, yn = true;
    bool zp = true, zn = true;

    for (uint8_t j = 0; j < 3; ++j) {
      points[j] = verts[indices[i*3+j]] - mOrigin;
      xp = xp && points[j].x > mResolution;
      xn = xn && points[j].x < 0;
      yp = yp && points[j].y > mResolution;
      yn = yn && points[j].y < 0;
      zp = zp && points[j].z > mResolution;
      zn = zn && points[j].z < 0;
    }

    // Check if all points are to one side
    if (xp || xn || yp || yn || zp || zn) continue;

    voxelizeTriangle(points);
  }
}
void VoxelGrid::writeToFile(const std::string& pPath) {
  std::ofstream fout;
  openFileWrite(fout, pPath);
  writeMetaData(fout);

  fout.write(&mVoxelData.at(0), mVoxelData.size());

  fout.close();
}

void VoxelGrid::writeToFileCompressed(const std::string& pPath, uint64_t* pVoxelsComplete) {
  std::ofstream fout;
  openFileWrite(fout, pPath);
  writeMetaData(fout);

  std::vector<uint64_t> compressedData = generateCompressedVoxelData(pVoxelsComplete);
  fout.write(reinterpret_cast<char*>(&compressedData.at(0)), compressedData.size() * sizeof(uint64_t));

  fout.close();
}

void VoxelGrid::loadFromFile(const std::string& pPath) {
  // Open file
  std::ifstream fin;
  fin.open(pPath, std::ios::binary | std::ios::in);
  
  // Read resolution
  fin.read(reinterpret_cast<char*>(&mResolution), sizeof(uint32_t));

  init();
  
  // Load
  // mVoxelData.resize(mVolume / 8);
  fin.read(&mVoxelData.at(0), mVolume / 8);
  // char byte;
  // std::println("Loaded file in {0}, decompressing...", t.getTime());
  // uint i = 0;
  // glm::uvec3 pos;
  // for (char byte: mVoxelData) {
  // // while (fin.read(&byte, 1)) {
  //   for (uint bitIndex = 0; bitIndex < 8; ++bitIndex, ++i, byte >>= 1) {
  //     if (byte & 1) {
  //       pos.z = i / (mResolution * mResolution);
  //       pos.y = (i / mResolution) % mResolution;
  //       pos.x = i % mResolution;
  //       mVoxelGrid[pos.x][pos.y][pos.z] = true;
  //     }
  //   }
  // }

  fin.close();
}

void VoxelGrid::loadFromFileCompressed(const std::string& pPath) {
  return;
}

void VoxelGrid::setLogStream(std::ostream* pStream, std::mutex* pMutex) {
  if (!pMutex)
    mLogMutex = &mDefaultLogMutex;
  else
    mLogMutex = pMutex;
  mLogStream = pStream;
}

int VoxelGrid::insert(const glm::uvec3& pPos, const insertFunc_t& pInsertFunc) {
  if (pPos.x < 0 || pPos.x >= mResolution ||
      pPos.y < 0 || pPos.y >= mResolution ||
      pPos.z < 0 || pPos.z >= mResolution) {
    std::lock_guard<std::mutex> lock(*mLogMutex);
    std::println(*mLogStream, "Couldn't insert voxel: {}, {}, {}. Reason: voxel outside grid", pPos.x, pPos.y, pPos.z);
    return 1;
  }

  uint64_t globalIndex = zorder(pPos);
  // uint globalIndex = pPos.x + pPos.y * mResolution + pPos.z * mResolution * mResolution;
  uint64_t byteIndex = globalIndex >> 3;
  uint localIndex = globalIndex & 0b111;
  char mask = 1 << localIndex;
  if (!(mVoxelData[byteIndex] & mask)) ++mVoxelCount;
  mVoxelData[byteIndex] |= mask;
  return 0;
}

void VoxelGrid::voxelizeTriangle(std::array<glm::vec3, 3> pPoints) {
  // Check if all points are the same
  if (glm::floor(pPoints[0]) == glm::floor(pPoints[1]) && glm::floor(pPoints[0]) == glm::floor(pPoints[2])) {
    insert(glm::floor(pPoints[0]));
    return;
  }

  // Find dominant axis
  glm::vec3 minPoint = glm::min(pPoints[0], glm::min(pPoints[1], pPoints[2]));
  glm::vec3 maxPoint = glm::max(pPoints[0], glm::max(pPoints[1], pPoints[2]));
  glm::vec3 boxSize = maxPoint - minPoint;

  uint8_t dominantAxisIndex = 0;
  if ((boxSize.x >= boxSize.y) && (boxSize.x >= boxSize.z))      dominantAxisIndex = 0;
  else if ((boxSize.y >= boxSize.x) && (boxSize.y >= boxSize.z)) dominantAxisIndex = 1;
  else if ((boxSize.z >= boxSize.x) && (boxSize.z >= boxSize.y)) dominantAxisIndex = 2;
  std::array<uint8_t, 2> nonDominantAxisIndices;
  nonDominantAxisIndices[0] = (dominantAxisIndex + 1) % 3;
  nonDominantAxisIndices[1] = (dominantAxisIndex + 2) % 3;

  // Sort by dominant axis
  if (v3index(pPoints[0], dominantAxisIndex) > v3index(pPoints[1], dominantAxisIndex)) std::swap(pPoints[0], pPoints[1]);
  if (v3index(pPoints[1], dominantAxisIndex) > v3index(pPoints[2], dominantAxisIndex)) std::swap(pPoints[1], pPoints[2]);
  if (v3index(pPoints[0], dominantAxisIndex) > v3index(pPoints[1], dominantAxisIndex)) std::swap(pPoints[0], pPoints[1]);

  // Calculate direction and inverse for lines
  glm::vec3 lhDir = pPoints[2] - pPoints[0];
  glm::vec3 lhDirInv = 1.f/lhDir;

  glm::vec3 lmDir = pPoints[1] - pPoints[0];
  glm::vec3 lmDirInv = 1.f/lmDir;

  glm::vec3 mhDir = pPoints[2] - pPoints[1];
  glm::vec3 mhDirInv = 1.f/mhDir;

  // Iterate over along the dominant axis traversing between the two lines
  for (float dominantAxisValue = std::max(0.f, v3index(pPoints[0], dominantAxisIndex)); dominantAxisValue <= std::min((float)mResolution, v3index(pPoints[2], dominantAxisIndex)); ++dominantAxisValue) {
    glm::vec3* l1Dir = &lmDir;
    glm::vec3* l1DirInv = &lmDirInv;
    glm::vec3* l2Dir = &lhDir;
    glm::vec3* l2DirInv = &lhDirInv;
    if ((v3index(pPoints[1], dominantAxisIndex) < dominantAxisValue) || (v3index(lmDir, dominantAxisIndex) == 0)) {
      l1Dir = &mhDir;
      l1DirInv = &mhDirInv;
    }

    glm::vec2 l1Pos;
    glm::vec2 l2Pos;

    // TODO: Use last save second intersection positions for the next iteration
    
    // L2
    // Normal plane
    float t = (dominantAxisValue - v3index(pPoints[0], dominantAxisIndex)) * v3index(*l2DirInv, dominantAxisIndex);
    glm::vec2 l2Posa = glm::vec2(v3index(pPoints[0], nonDominantAxisIndices[0]), v3index(pPoints[0], nonDominantAxisIndices[1])) // Origin
                      + t * glm::vec2(v3index(*l2Dir, nonDominantAxisIndices[0]), v3index(*l2Dir, nonDominantAxisIndices[1]));
    // Next plane
    t = (std::min(dominantAxisValue + 1, v3index(pPoints[2], dominantAxisIndex)) - v3index(pPoints[0], dominantAxisIndex)) * v3index(*l2DirInv, dominantAxisIndex);
    glm::vec2 l2Posb = glm::vec2(v3index(pPoints[0], nonDominantAxisIndices[0]), v3index(pPoints[0], nonDominantAxisIndices[1])) // Origin
                      + t * glm::vec2(v3index(*l2Dir, nonDominantAxisIndices[0]), v3index(*l2Dir, nonDominantAxisIndices[1]));

    if ((v3index(pPoints[1], dominantAxisIndex) >= dominantAxisValue) && (v3index(pPoints[1], dominantAxisIndex) <= dominantAxisValue + 1)) {
      l1Pos = glm::vec2(v3index(pPoints[1], nonDominantAxisIndices[0]), v3index(pPoints[1], nonDominantAxisIndices[1]));

      // Choose which has the greatest length
      if (glm::length(l2Posa - l1Pos) >= glm::length(l2Posb - l1Pos)) l2Pos = l2Posa;
      else                                                            l2Pos = l2Posb;
    }
    else {
      // L1
      // Normal plane
      t = (dominantAxisValue - v3index(pPoints[1], dominantAxisIndex)) * v3index(*l1DirInv, dominantAxisIndex);
      glm::vec2 l1Posa = glm::vec2(v3index(pPoints[1], nonDominantAxisIndices[0]), v3index(pPoints[1], nonDominantAxisIndices[1])) // Origin
                        + t * glm::vec2(v3index(*l1Dir, nonDominantAxisIndices[0]), v3index(*l1Dir, nonDominantAxisIndices[1]));
      // Next plane
      t = (std::min(dominantAxisValue + 1, v3index(pPoints[2], dominantAxisIndex)) - v3index(pPoints[1], dominantAxisIndex)) * v3index(*l1DirInv, dominantAxisIndex);
      glm::vec2 l1Posb = glm::vec2(v3index(pPoints[1], nonDominantAxisIndices[0]), v3index(pPoints[1], nonDominantAxisIndices[1])) // Origin
                        + t * glm::vec2(v3index(*l1Dir, nonDominantAxisIndices[0]), v3index(*l1Dir, nonDominantAxisIndices[1]));

      // Choose which has the greatest length
      if (glm::length(l2Posa - l1Posa) >= glm::length(l2Posb - l1Posa)) l2Pos = l2Posa;
      else                                                              l2Pos = l2Posb;

      if (glm::length(l2Pos - l1Posa) >= glm::length(l2Pos - l1Posb)) l1Pos = l1Posa;
      else                                                            l1Pos = l1Posb;
    }

    // TODO: This could be optimized by mixing it in rahter than just doing it after
    uint index = 4;
    if ((v3index(pPoints[0], dominantAxisIndex) >= dominantAxisValue) && (v3index(pPoints[0], dominantAxisIndex) <= dominantAxisValue + 1))
      index = 0;
    if ((v3index(pPoints[2], dominantAxisIndex) >= dominantAxisValue) && (v3index(pPoints[2], dominantAxisIndex) <= dominantAxisValue + 1))
      index = 2;
    if (index != 4) {
      glm::vec2 v2(v3index(pPoints[index], nonDominantAxisIndices[0]), v3index(pPoints[index], nonDominantAxisIndices[1]));
      if (glm::length(l2Pos - v2) >= glm::length(l2Pos - l1Pos))      l1Pos = v2;
      else if (glm::length(l1Pos - v2) >= glm::length(l1Pos - l2Pos)) l2Pos = v2;
    }

    glm::vec2 dir = l2Pos - l1Pos; // Dir from l1 to l2 intersection pPoints
    glm::vec2 dirInv = 1.f/dir;

    drawLine2(dominantAxisIndex, dominantAxisValue, l1Pos, l2Pos, dir, dirInv);
  }
}

uint64_t VoxelGrid::getVoxelCount() {
  return mVoxelCount;
}

uint64_t VoxelGrid::getVolume() {
  return mVolume;
}

uint VoxelGrid::getResolution() {
  return mResolution;
}

float VoxelGrid::getMaxDepth() {
  return mMaxDepth;
}

char* VoxelGrid::getVoxelDataByte(const glm::uvec3& pPos) {
  uint64_t globalIndex = zorder(pPos);
  uint64_t byteIndex = globalIndex >> 3;
  uint localIndex = globalIndex & 0b111;
  return &mVoxelData[byteIndex];
}

const bool VoxelGrid::queryVoxel(const glm::uvec3& pPos) {
  if (pPos.x < 0 || pPos.x >= mResolution ||
      pPos.y < 0 || pPos.y >= mResolution ||
      pPos.z < 0 || pPos.z >= mResolution)
    throw 1;
  
  uint64_t globalIndex = zorder(pPos);
  // uint globalIndex = pPos.x + pPos.y * mResolution + pPos.z * mResolution * mResolution;
  uint64_t byteIndex = globalIndex >> 3;
  uint localIndex = globalIndex & 0b111;
  char mask = 1 << localIndex;
  return mVoxelData[byteIndex] & mask;
}

// const std::vector<std::vector<std::vector<bool>>>& VoxelGrid::getVoxelData() {
  // return mVoxelGrid;
// }

const std::vector<char>& VoxelGrid::getVoxelDataBits() {
  return mVoxelData;
}

std::vector<uint64_t> VoxelGrid::generateCompressedVoxelData(uint64_t* pVoxelsComplete) {
  std::vector<uint64_t> counts;
  uint64_t count = 0;
  bool value = false;
  for (uint z = 0; z < mResolution; ++z) {
    for (uint y = 0; y < mResolution; ++y) {
      for (uint x = 0; x < mResolution; ++x) {
        if (queryVoxel(glm::vec3(x, y, z)) != value) {
          // Write count
          counts.push_back(count);
          value = !value;
          count = 0;
        }
        ++count;
        ++*pVoxelsComplete;
      }
    }
  }
  return counts;
}

void VoxelGrid::init() {
  uint64_t res = mResolution;
  mVolume = res * res * res;
  mMaxDepth = std::log2f(mResolution);

  // mVoxelGrid.clear();
  mVoxelData.clear();

  // Init mVoxelGrid and mVoxelData
  // mVoxelGrid = std::vector<std::vector<std::vector<bool>>>(mResolution, std::vector<std::vector<bool>>(mResolution, std::vector<bool>(mResolution, false)));

  mVoxelData = std::vector<char>(mVolume >> 3, 0);
}

void VoxelGrid::openFileWrite(std::ofstream& pFout, const std::string& pPath) {
  pFout.open(pPath, std::ios::out | std::ios::binary);
  if (!pFout.is_open())
    throw "Could not open output file";
}

void VoxelGrid::writeMetaData(std::ofstream& pFout) {
  // Write resolution
  pFout.write(reinterpret_cast<char*>(&mResolution), sizeof(uint32_t));
}

void VoxelGrid::drawLine2(uint8_t pDominantAxisIndex, float pDominantAxisValue, const glm::vec2& pStart, const glm::vec2& pEnd, const glm::vec2& pDir, const glm::vec2& pDirInv) {
  // std::println("Started drawLine2");
  glm::vec3 v = glm::floor(toVec3(pDominantAxisValue, pStart.x, pStart.y, pDominantAxisIndex));

  insert(v);

  glm::vec2 voxelPos = glm::floor(pStart);

  while (  voxelPos != glm::floor(pEnd) &&
         !(voxelPos.x > mResolution && pDir.x > 0) &&
         !(voxelPos.x < 0 && pDir.x < 0) &&
         !(voxelPos.y > mResolution && pDir.y > 0) &&
         !(voxelPos.y < 0 && pDir.y < 0)) {
    // std::println("voxelPos: ({}, {}), l1Pos: ({}, {}), l2pos: ({}, {}), signDirInv: ({}, {})", voxelPos.x, voxelPos.y, pStart.x, pStart.y, pEnd.x, pEnd.y, glm::sign(pDirInv.x), glm::sign(pDirInv.y));
    float plane1 = voxelPos.x + std::max(0.f, glm::sign(pDirInv.x));
    float plane2 = voxelPos.y + std::max(0.f, glm::sign(pDirInv.y));

    if (plane1 == pEnd.x && plane2 == pEnd.y) { // If the destination point lies exactly on the planes we increment both x and y
      voxelPos += glm::sign(pDirInv);
      v = glm::floor(toVec3(pDominantAxisValue, voxelPos.x, voxelPos.y, pDominantAxisIndex));
      insert(v);
      break;
    }
    
    float t1 = (plane1 - pStart.x) * pDirInv.x;
    float t2 = (plane2 - pStart.y) * pDirInv.y;

    if (t1 < t2) voxelPos.x += glm::sign(pDirInv.x);
    else         voxelPos.y += glm::sign(pDirInv.y);

    v = glm::floor(toVec3(pDominantAxisValue, voxelPos.x, voxelPos.y, pDominantAxisIndex));

    insert(v);
  }

  // std::println("Started");
  // int steps;
  // steps = std::max(abs(pDir.x), abs(pDir.y));
  // std::println("Steps: {}, dir: {}, {}, ceil(dir): {}, {}", steps, pDir.x, pDir.y, abs(pDir.x), abs(pDir.y));
  // glm::vec2 incr = pDir / (float)steps;
  // glm::vec2 voxelPos = pStart;
  // for (uint i = 0; i < steps; ++i) {
  //   glm::vec3 v = toVec3(pDominantAxisValue, glm::floor(voxelPos.x), glm::floor(voxelPos.y), pDominantAxisIndex);
  //   // glm::vec3 v = toVec3(a, floor(voxelPos.x), floor(voxelPos.y), dominantAxisIndex);
  //   if (!mVoxelGrid[v.x][v.y][v.z]) {
  //     std::println("Inserting");
  //     insert(v, {glm::vec4(redValue, 0, pBlueValue, 0)});
  //     redValue += 0.1f;
  //     if (redValue > 1.f)
  //       redValue = 0.f;
  //     mVoxelGrid[v.x][v.y][v.z] = true;
  //   }
  //   voxelPos += incr;
  // }
  // std::println("pEnd: {}, {}", pEnd.x, pEnd.y);
  // glm::vec3 v = glm::floor(toVec3(pDominantAxisValue, pEnd.x, pEnd.y, pDominantAxisIndex));
  // if (!mVoxelGrid[v.x][v.y][v.z]) {
  //   ++mVoxelCount;
  //   insert(v, {glm::vec4(0, 0, 1, 0)});
  //   mVoxelGrid[v.x][v.y][v.z] = true;
  // }
  // std::println("Finished");
}

glm::vec3 VoxelGrid::toVec3(float a, float b, float c, uint8_t pDominantAxisIndex) {
  glm::vec3 v;
  if (pDominantAxisIndex == 0) {
    v.x = a;
    v.y = b;
    v.z = c;
  }
  else if (pDominantAxisIndex == 1) {
    v.x = c;
    v.y = a;
    v.z = b;
  }
  else if (pDominantAxisIndex == 2) {
    v.x = b;
    v.y = c;
    v.z = a;
  }
  return v;
}

uint64_t VoxelGrid::zorder(const glm::uvec3& pPos) {
  uint64_t x = pPos.x, y = pPos.y, z = pPos.z;
  static const uint64_t B[] = {0x00000000FF0000FF, 0x000000F00F00F00F,
                               0x00000C30C30C30C3, 0X0000249249249249};
  static const int S[] =  {16, 8, 4, 2}; 
  static const uint64_t MAXINPUT = 65536;

  assert( (x < MAXINPUT) && (y < MAXINPUT) && (z < MAXINPUT) );

  x = (x | (x << S[0])) & B[0];
  x = (x | (x << S[1])) & B[1];
  x = (x | (x << S[2])) & B[2];
  x = (x | (x << S[3])) & B[3];

  y = (y | (y << S[0])) & B[0];
  y = (y | (y << S[1])) & B[1];
  y = (y | (y << S[2])) & B[2];
  y = (y | (y << S[3])) & B[3];

  z = (z | (z <<  S[0])) & B[0];
  z = (z | (z <<  S[1])) & B[1];
  z = (z | (z <<  S[2])) & B[2];
  z = (z | (z <<  S[3])) & B[3];

  uint64_t result = x | (y << 1) | (z << 2);
  return result;
}

