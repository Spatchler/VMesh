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

void VoxelGrid::clear() {
  mVoxelCount = 0;
  std::fill(mVoxelData.begin(), mVoxelData.end(), 0);
}

void VoxelGrid::IntersectVoxelizeMesh(Mesh& pMesh, uint* pTrisComplete) {
  const std::vector<Vertex>& verts = pMesh.getVertices();
  const std::vector<uint>& indices = pMesh.getIndices();

  std::array<glm::vec3, 3> points;

  for (uint i = 0; i < pMesh.getTriCount(); ++i, ++*pTrisComplete) { // Every triangle
    // Add points to array
    for (uint j = 0; j < 3; ++j) points[j] = verts[indices[i*3+j]].pos - mOrigin;

    glm::vec3 minPoint = glm::min(points[0], glm::min(points[1], points[2]));
    glm::vec3 maxPoint = glm::max(points[0], glm::max(points[1], points[2]));

    minPoint = glm::floor(minPoint);
    maxPoint = glm::ceil(maxPoint);

    minPoint = glm::max(minPoint, glm::vec3(0));
    maxPoint = glm::min(maxPoint, glm::vec3(mResolution));

    for (glm::uvec3 pos = minPoint; pos.z <= maxPoint.z; ++pos.z) for (pos.y = minPoint.y; pos.y <= maxPoint.y; ++pos.y) for (pos.x = minPoint.x; pos.x <= maxPoint.x; ++pos.x) // Every block
      // if (!queryVoxel(pos)) if (triangleBoxIntersect(points, glm::vec3(pos) + 0.5f)) insert(pos);
      if (triangleBoxIntersect(points, glm::vec3(pos) + 0.5f)) insert(pos, 1);
  }
}

void VoxelGrid::IntersectVoxelizeModel(Model& pModel, uint* pTrisComplete) {
  for (uint i = 0; i < pModel.getNumMeshes(); ++i)
    IntersectVoxelizeMesh(pModel.getMesh(i), pTrisComplete);
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
    v = pVert[i];
    if (pNormal[i] > 0.0f) {
      vmin[i] = -0.5f - v;
      vmax[i] =  0.5f - v;
    }
    else {
      vmin[i] =  0.5f - v;
      vmax[i] = -0.5f - v;
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

void VoxelGrid::DDAvoxelizeMesh(Mesh& pMesh, uint* pTrisComplete, Texture* pTex) {
  const std::vector<Vertex>& verts = pMesh.getVertices();
  const std::vector<uint>& indices = pMesh.getIndices();

  for (uint i = 0; i < pMesh.getTriCount(); ++i, ++*pTrisComplete) { // Every triangle
    std::array<Vertex, 3> points;

    bool xp = true, xn = true;
    bool yp = true, yn = true;
    bool zp = true, zn = true;

    for (uint8_t j = 0; j < 3; ++j) {
      points[j] = verts[indices[i*3+j]];
      points[j].pos -= mOrigin;
      xp = xp && points[j].pos.x > mResolution;
      xn = xn && points[j].pos.x < 0;
      yp = yp && points[j].pos.y > mResolution;
      yn = yn && points[j].pos.y < 0;
      zp = zp && points[j].pos.z > mResolution;
      zn = zn && points[j].pos.z < 0;
    }

    // Check if all points are to one side
    if (xp || xn || yp || yn || zp || zn) continue;

    DDAvoxelizeTriangle(points, pTex);
  }
}

void VoxelGrid::DDAvoxelizeModel(Model& pModel, uint* pTrisComplete, bool pColoured) {
  for (uint i = 0; i < pModel.getNumMeshes(); ++i) {
    Mesh& m = pModel.getMesh(i);
    Texture* t = NULL;
    if (pColoured) t = &pModel.getDiffuseMap(m.getMatIndex());
    DDAvoxelizeMesh(m, pTrisComplete, t);
  }
}

void VoxelGrid::writeToFile(const std::string& pPath) {
  std::ofstream fout;
  fout.open(pPath, std::ios::out | std::ios::binary);
  if (!fout.is_open()) throw std::runtime_error("Could not open output file");

  writeMetaData(fout);

  fout.write(reinterpret_cast<char*>(&mVoxelData.at(0)), mVoxelData.size());

  fout.close();
}

void VoxelGrid::writeToFileCompressed(const std::string& pPath, uint64_t* pVoxelsComplete) {
  std::ofstream fout;
  fout.open(pPath, std::ios::out | std::ios::binary);
  if (!fout.is_open()) throw std::runtime_error("Could not open output file");

  writeMetaData(fout);

  std::vector<uint64_t> compressedData = generateCompressedVoxelData(pVoxelsComplete);

  uint32_t countsLen = compressedData.size();
  fout.write(reinterpret_cast<char*>(&countsLen), sizeof(uint32_t));

  fout.write(reinterpret_cast<char*>(&compressedData.at(0)), compressedData.size() * sizeof(uint64_t));

  fout.close();
}

void VoxelGrid::loadFromFile(const std::string& pPath) {
  // Open file
  std::ifstream fin;
  fin.open(pPath, std::ios::binary | std::ios::in);
  if (!fin.is_open()) throw std::runtime_error("Could not open input file");

  readMetaData(fin);

  init();
  
  // Load
  fin.read(reinterpret_cast<char*>(&mVoxelData.at(0)), mVolume);

  fin.close();
}

void VoxelGrid::loadFromFileCompressed(const std::string& pPath) {
  // Open file
  std::ifstream fin;
  fin.open(pPath, std::ios::binary | std::ios::in);
  if (!fin.is_open()) throw std::runtime_error("Could not open input file");

  readMetaData(fin);

  // Read counts len
  uint32_t countsLen;
  fin.read(reinterpret_cast<char*>(&countsLen), sizeof(uint32_t));

  init();

  // Load
  std::vector<uint64_t> counts;
  counts.resize(countsLen);
  fin.read(reinterpret_cast<char*>(&counts.at(0)), countsLen * sizeof(uint64_t));

  bool isAir = true;

  uint64_t index = 0;
  for (uint64_t i: counts) {
    if (!isAir) for (uint64_t j = 0; j < i; ++j) insert(index + j, 1);
    index += i;
    isAir = !isAir;
  }

  fin.close();
}

void VoxelGrid::setLogStream(std::ostream* pStream, std::mutex* pMutex) {
  if (!pMutex) mLogMutex = &mDefaultLogMutex;
  else         mLogMutex = pMutex;
  mLogStream = pStream;
}

void VoxelGrid::insert(uint64_t pIndex, uint8_t pCol) {
  if (mVoxelData.at(pIndex) == 0) ++mVoxelCount;
  mVoxelData.at(pIndex) = pCol;
  if (pCol == 0) --mVoxelCount;
}

void VoxelGrid::insert(const glm::uvec3& pPos, uint8_t pCol) {
  insert(zorder(pPos), pCol);
}

void VoxelGrid::DDAvoxelizeTriangle(std::array<Vertex, 3> pVerts, Texture* pTex) {
  uint8_t col = 1;
  if (pTex) col = mPalette.addColour(pTex->sample(pVerts[0].texCoord), 0.01);
  // Check if all points are the same
  if (glm::floor(pVerts[0].pos) == glm::floor(pVerts[1].pos) && glm::floor(pVerts[0].pos) == glm::floor(pVerts[2].pos)) {
    insert(glm::floor(pVerts[0].pos), col);
    return;
  }

  // Find dominant axis
  glm::vec3 minPoint = glm::min(pVerts[0].pos, glm::min(pVerts[1].pos, pVerts[2].pos));
  glm::vec3 maxPoint = glm::max(pVerts[0].pos, glm::max(pVerts[1].pos, pVerts[2].pos));
  glm::vec3 boxSize = maxPoint - minPoint;

  uint8_t dominantAxisIndex = 0;
  if ((boxSize.x >= boxSize.y) && (boxSize.x >= boxSize.z))      dominantAxisIndex = 0;
  else if ((boxSize.y >= boxSize.x) && (boxSize.y >= boxSize.z)) dominantAxisIndex = 1;
  else if ((boxSize.z >= boxSize.x) && (boxSize.z >= boxSize.y)) dominantAxisIndex = 2;
  std::array<uint8_t, 2> nonDominantAxisIndices;
  nonDominantAxisIndices[0] = (dominantAxisIndex + 1) % 3;
  nonDominantAxisIndices[1] = (dominantAxisIndex + 2) % 3;

  // Sort by dominant axis
  if (pVerts[0].pos[dominantAxisIndex] > pVerts[1].pos[dominantAxisIndex]) std::swap(pVerts[0].pos, pVerts[1].pos);
  if (pVerts[1].pos[dominantAxisIndex] > pVerts[2].pos[dominantAxisIndex]) std::swap(pVerts[1].pos, pVerts[2].pos);
  if (pVerts[0].pos[dominantAxisIndex] > pVerts[1].pos[dominantAxisIndex]) std::swap(pVerts[0].pos, pVerts[1].pos);

  // Calculate direction and inverse for lines
  glm::vec3 lhDir = pVerts[2].pos - pVerts[0].pos;
  glm::vec3 lhDirInv = 1.f/lhDir;

  glm::vec3 lmDir = pVerts[1].pos - pVerts[0].pos;
  glm::vec3 lmDirInv = 1.f/lmDir;

  glm::vec3 mhDir = pVerts[2].pos - pVerts[1].pos;
  glm::vec3 mhDirInv = 1.f/mhDir;

  // Iterate over along the dominant axis traversing between the two lines
  for (float dominantAxisValue = std::max(0.f, pVerts[0].pos[dominantAxisIndex]); dominantAxisValue <= std::min((float)mResolution, pVerts[2].pos[dominantAxisIndex]); ++dominantAxisValue) {
    glm::vec3* l1Dir = &lmDir;
    glm::vec3* l1DirInv = &lmDirInv;
    glm::vec3* l2Dir = &lhDir;
    glm::vec3* l2DirInv = &lhDirInv;
    if ((pVerts[1].pos[dominantAxisIndex] < dominantAxisValue) || (lmDir[dominantAxisIndex] == 0)) {
      l1Dir = &mhDir;
      l1DirInv = &mhDirInv;
    }

    glm::vec2 l1Pos;
    glm::vec2 l2Pos;

    // TODO: Use last save second intersection positions for the next iteration
    
    // L2
    // Normal plane
    float t = (dominantAxisValue - pVerts[0].pos[dominantAxisIndex]) * (*l2DirInv)[dominantAxisIndex];
    glm::vec2 l2Posa = glm::vec2(pVerts[0].pos[nonDominantAxisIndices[0]], pVerts[0].pos[nonDominantAxisIndices[1]]) // Origin
                      + t * glm::vec2((*l2Dir)[nonDominantAxisIndices[0]], (*l2Dir)[nonDominantAxisIndices[1]]);
    // Next plane
    t = (std::min(dominantAxisValue + 1, pVerts[2].pos[dominantAxisIndex]) - pVerts[0].pos[dominantAxisIndex]) * (*l2DirInv)[dominantAxisIndex];
    glm::vec2 l2Posb = glm::vec2(pVerts[0].pos[nonDominantAxisIndices[0]], pVerts[0].pos[nonDominantAxisIndices[1]]) // Origin
                      + t * glm::vec2((*l2Dir)[nonDominantAxisIndices[0]], (*l2Dir)[nonDominantAxisIndices[1]]);

    if ((pVerts[1].pos[dominantAxisIndex] >= dominantAxisValue) && (pVerts[1].pos[dominantAxisIndex] <= dominantAxisValue + 1)) {
      l1Pos = glm::vec2(pVerts[1].pos[nonDominantAxisIndices[0]], pVerts[1].pos[nonDominantAxisIndices[1]]);

      // Choose which has the greatest length
      if (glm::length(l2Posa - l1Pos) >= glm::length(l2Posb - l1Pos)) l2Pos = l2Posa;
      else                                                            l2Pos = l2Posb;
    }
    else {
      // L1
      // Normal plane
      t = (dominantAxisValue - pVerts[1].pos[dominantAxisIndex]) * (*l1DirInv)[dominantAxisIndex];
      glm::vec2 l1Posa = glm::vec2(pVerts[1].pos[nonDominantAxisIndices[0]], pVerts[1].pos[nonDominantAxisIndices[1]]) // Origin
                        + t * glm::vec2((*l1Dir)[nonDominantAxisIndices[0]], (*l1Dir)[nonDominantAxisIndices[1]]);
      // Next plane
      t = (std::min(dominantAxisValue + 1, pVerts[2].pos[dominantAxisIndex]) - pVerts[1].pos[dominantAxisIndex]) * (*l1DirInv)[dominantAxisIndex];
      glm::vec2 l1Posb = glm::vec2(pVerts[1].pos[nonDominantAxisIndices[0]], pVerts[1].pos[nonDominantAxisIndices[1]]) // Origin
                        + t * glm::vec2((*l1Dir)[nonDominantAxisIndices[0]], (*l1Dir)[nonDominantAxisIndices[1]]);

      // Choose which has the greatest length
      if (glm::length(l2Posa - l1Posa) >= glm::length(l2Posb - l1Posa)) l2Pos = l2Posa;
      else                                                              l2Pos = l2Posb;

      if (glm::length(l2Pos - l1Posa) >= glm::length(l2Pos - l1Posb)) l1Pos = l1Posa;
      else                                                            l1Pos = l1Posb;
    }

    // TODO: This could be optimized by mixing it in rahter than just doing it after
    uint index = 4;
    if ((pVerts[0].pos[dominantAxisIndex] >= dominantAxisValue) && (pVerts[0].pos[dominantAxisIndex] <= dominantAxisValue + 1))
      index = 0;
    if ((pVerts[2].pos[dominantAxisIndex] >= dominantAxisValue) && (pVerts[2].pos[dominantAxisIndex] <= dominantAxisValue + 1))
      index = 2;
    if (index != 4) {
      glm::vec2 v2(pVerts[index].pos[nonDominantAxisIndices[0]], pVerts[index].pos[nonDominantAxisIndices[1]]);
      if (glm::length(l2Pos - v2) >= glm::length(l2Pos - l1Pos))      l1Pos = v2;
      else if (glm::length(l1Pos - v2) >= glm::length(l1Pos - l2Pos)) l2Pos = v2;
    }

    glm::vec2 dir = l2Pos - l1Pos; // Dir from l1 to l2 intersection pPoints
    glm::vec2 dirInv = 1.f/dir;
    
    drawLine2(dominantAxisIndex, dominantAxisValue, l1Pos, l2Pos, dir, dirInv, pTex, col);
  }
}

bool VoxelGrid::isRegionAllSame(const glm::uvec3& pOrigin, uint pSize) {
  uint8_t first = queryVoxelData(pOrigin);
  for (glm::uvec3 p(0); p.z < pSize; ++p.z) for (p.y = 0; p.y < pSize; ++p.y) for (p.x = 0; p.x < pSize; ++p.x)
    if (queryVoxelData(pOrigin + p) != first) return false;
  return true;
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

uint8_t VoxelGrid::queryVoxelData(uint64_t pIndex) {
  return mVoxelData.at(pIndex);
}

uint8_t VoxelGrid::queryVoxelData(const glm::uvec3& pPos) {
  return queryVoxelData(zorder(pPos));
}

const std::vector<uint8_t>& VoxelGrid::getVoxelData() {
  return mVoxelData;
}

std::vector<uint64_t> VoxelGrid::generateCompressedVoxelData(uint64_t* pVoxelsComplete) {
  std::vector<uint64_t> counts;
  uint64_t count = 0;
  uint8_t value = 0;
  for (uint8_t i: mVoxelData) {
    if (i == value) {
      ++count;
      continue;
    }
    counts.push_back(value);
    counts.push_back(count);
    value = i;
    count = 1;
  }

  if (value) {
    counts.push_back(value);
    counts.push_back(count);
  }

  return counts;
}

void VoxelGrid::init() {
  uint64_t res = mResolution;
  mResolution2 = res * res;
  mVolume = mResolution2 * res;
  mMaxDepth = std::log2f(mResolution);

  mVoxelData.resize(mVolume, 0);
}

void VoxelGrid::writeMetaData(std::ofstream& pFout) {
  pFout.write(reinterpret_cast<char*>(&mResolution), sizeof(uint32_t));
  pFout.write(reinterpret_cast<char*>(&mVoxelCount), sizeof(uint64_t));
}

void VoxelGrid::readMetaData(std::ifstream& pFin) {
  pFin.read(reinterpret_cast<char*>(&mResolution), sizeof(uint32_t));
  pFin.read(reinterpret_cast<char*>(&mVoxelCount), sizeof(uint64_t));
}

void VoxelGrid::drawLine2(uint8_t pDominantAxisIndex, float pDominantAxisValue, const glm::vec2& pStart, const glm::vec2& pEnd, const glm::vec2& pDir, const glm::vec2& pDirInv, Texture* pTex, uint8_t pCol) {
  // std::println("Started drawLine2");
  glm::vec3 v = glm::floor(toVec3(pDominantAxisValue, pStart.x, pStart.y, pDominantAxisIndex));

  insert(v, pCol);

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
      insert(v, pCol);
      break;
    }
    
    float t1 = (plane1 - pStart.x) * pDirInv.x;
    float t2 = (plane2 - pStart.y) * pDirInv.y;

    if (t1 < t2) voxelPos.x += glm::sign(pDirInv.x);
    else         voxelPos.y += glm::sign(pDirInv.y);

    v = glm::floor(toVec3(pDominantAxisValue, voxelPos.x, voxelPos.y, pDominantAxisIndex));

    insert(v, pCol);
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
  v[pDominantAxisIndex] = a;
  v[(pDominantAxisIndex + 1) % 3] = b;
  v[(pDominantAxisIndex + 2) % 3] = c;
  // glm::vec3 v;
  // if (pDominantAxisIndex == 0) {
  //   v.x = a;
  //   v.y = b;
  //   v.z = c;
  // }
  // else if (pDominantAxisIndex == 1) {
  //   v.x = c;
  //   v.y = a;
  //   v.z = b;
  // }
  // else if (pDominantAxisIndex == 2) {
  //   v.x = b;
  //   v.y = c;
  //   v.z = a;
  // }
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

