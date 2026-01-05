#include "voxelGrid.hpp"

using namespace VMesh;

VoxelGrid::VoxelGrid(uint pResolution)
:mResolution(pResolution) {
  mLogStream = &mDefaultLogStream;
  mLogMutex = &mDefaultLogMutex;
  init();
}

void VoxelGrid::voxelizeMesh(Mesh& pMesh, uint* pTrisComplete, const insertFunc_t& pInsertFunc) {
  const std::vector<glm::vec3>& verts = pMesh.getVertices();
  const std::vector<uint>& indices = pMesh.getIndices();

  for (uint i = 0; i < pMesh.getTriCount(); ++i) { // Every triangle
    // Add points to array
    std::array<glm::vec3, 3> points;
    for (uint j = 0; j < 3; ++j)
      points[j] = verts[indices[i*3+j]];

    // Find dominant axis
    glm::vec3 minPoint = glm::min(points[0], glm::min(points[1], points[2]));
    glm::vec3 maxPoint = glm::max(points[0], glm::max(points[1], points[2]));
    glm::vec3 boxSize = maxPoint - minPoint;

    uint8_t dominantAxisIndex = 0;
    if ((boxSize.x >= boxSize.y) && (boxSize.x >= boxSize.z))
      dominantAxisIndex = 0;
    else if ((boxSize.y >= boxSize.x) && (boxSize.y >= boxSize.z))
      dominantAxisIndex = 1;
    else if ((boxSize.z >= boxSize.x) && (boxSize.z >= boxSize.y))
      dominantAxisIndex = 2;
    std::array<uint8_t, 2> nonDominantAxisIndices;
    nonDominantAxisIndices[0] = (dominantAxisIndex + 1) % 3;
    nonDominantAxisIndices[1] = (dominantAxisIndex + 2) % 3;

    // glm::vec3 normal = glm::cross(points[1] - points[0], points[2] - points[0]);
    // normal = glm::normalize(normal);

    // Sort by dominant axis
    if (v3index(points[0], dominantAxisIndex) > v3index(points[1], dominantAxisIndex)) std::swap(points[0], points[1]);
    if (v3index(points[1], dominantAxisIndex) > v3index(points[2], dominantAxisIndex)) std::swap(points[1], points[2]);
    if (v3index(points[0], dominantAxisIndex) > v3index(points[1], dominantAxisIndex)) std::swap(points[0], points[1]);

    // Calculate direction and inverse for lines
    glm::vec3 lhDir = points[2] - points[0];
    glm::vec3 lhDirInv = 1.f/lhDir;

    glm::vec3 lmDir = points[1] - points[0];
    glm::vec3 lmDirInv = 1.f/lmDir;

    glm::vec3 mhDir = points[2] - points[1];
    glm::vec3 mhDirInv = 1.f/mhDir;

    // Iterate over along the dominant axis traversing between the two lines
    for (float dominantAxisValue = v3index(points[0], dominantAxisIndex); dominantAxisValue <= v3index(points[2], dominantAxisIndex); ++dominantAxisValue) {
      glm::vec3* l1Dir = &lmDir;
      glm::vec3* l1DirInv = &lmDirInv;
      glm::vec3* l2Dir = &lhDir;
      glm::vec3* l2DirInv = &lhDirInv;
      if ((v3index(points[1], dominantAxisIndex) < dominantAxisValue) || (v3index(lmDir, dominantAxisIndex) == 0)) {
        l1Dir = &mhDir;
        l1DirInv = &mhDirInv;
      }

      glm::vec2 l1Pos;
      glm::vec2 l2Pos;

      // TODO: Use last save second intersection positions for the next iteration

      // L2
      // Normal plane
      float t = (dominantAxisValue - v3index(points[0], dominantAxisIndex)) * v3index(*l2DirInv, dominantAxisIndex);
      glm::vec2 l2Posa = glm::vec2(v3index(points[0], nonDominantAxisIndices[0]), v3index(points[0], nonDominantAxisIndices[1])) // Origin
                        + t * glm::vec2(v3index(*l2Dir, nonDominantAxisIndices[0]), v3index(*l2Dir, nonDominantAxisIndices[1]));
      // Next plane
      t = (std::min(dominantAxisValue + 1, v3index(points[2], dominantAxisIndex)) - v3index(points[0], dominantAxisIndex)) * v3index(*l2DirInv, dominantAxisIndex);
      glm::vec2 l2Posb = glm::vec2(v3index(points[0], nonDominantAxisIndices[0]), v3index(points[0], nonDominantAxisIndices[1])) // Origin
                        + t * glm::vec2(v3index(*l2Dir, nonDominantAxisIndices[0]), v3index(*l2Dir, nonDominantAxisIndices[1]));

      if ((v3index(points[1], dominantAxisIndex) >= dominantAxisValue) && (v3index(points[1], dominantAxisIndex) <= dominantAxisValue + 1)) {
        l1Pos = glm::vec2(v3index(points[1], nonDominantAxisIndices[0]), v3index(points[1], nonDominantAxisIndices[1]));

        // Choose which has the greatest length
        if (glm::length(l2Posa - l1Pos) >= glm::length(l2Posb - l1Pos))
          l2Pos = l2Posa;
        else
          l2Pos = l2Posb;
      }
      else {
        // L1
        // Normal plane
        t = (dominantAxisValue - v3index(points[1], dominantAxisIndex)) * v3index(*l1DirInv, dominantAxisIndex);
        glm::vec2 l1Posa = glm::vec2(v3index(points[1], nonDominantAxisIndices[0]), v3index(points[1], nonDominantAxisIndices[1])) // Origin
                          + t * glm::vec2(v3index(*l1Dir, nonDominantAxisIndices[0]), v3index(*l1Dir, nonDominantAxisIndices[1]));
        // Next plane
        t = (std::min(dominantAxisValue + 1, v3index(points[2], dominantAxisIndex)) - v3index(points[1], dominantAxisIndex)) * v3index(*l1DirInv, dominantAxisIndex);
        glm::vec2 l1Posb = glm::vec2(v3index(points[1], nonDominantAxisIndices[0]), v3index(points[1], nonDominantAxisIndices[1])) // Origin
                          + t * glm::vec2(v3index(*l1Dir, nonDominantAxisIndices[0]), v3index(*l1Dir, nonDominantAxisIndices[1]));

        // Choose which has the greatest length
        if (glm::length(l2Posa - l1Posa) >= glm::length(l2Posb - l1Posa))
          l2Pos = l2Posa;
        else
          l2Pos = l2Posb;

        if (glm::length(l2Pos - l1Posa) >= glm::length(l2Pos - l1Posb))
          l1Pos = l1Posa;
        else
          l1Pos = l1Posb;
      }

      // TODO: This could be optimized by mixing it in rahter than just doing it after
      uint index = 4;
      if ((v3index(points[0], dominantAxisIndex) >= dominantAxisValue) && (v3index(points[0], dominantAxisIndex) <= dominantAxisValue + 1))
        index = 0;
      if ((v3index(points[2], dominantAxisIndex) >= dominantAxisValue) && (v3index(points[2], dominantAxisIndex) <= dominantAxisValue + 1))
        index = 2;
      if (index != 4) {
        glm::vec2 v2(v3index(points[index], nonDominantAxisIndices[0]), v3index(points[index], nonDominantAxisIndices[1]));
        if (glm::length(l2Pos - v2) >= glm::length(l2Pos - l1Pos))
          l1Pos = v2;
        else if (glm::length(l1Pos - v2) >= glm::length(l1Pos - l2Pos))
          l2Pos = v2;
      }

      glm::vec2 dir = l2Pos - l1Pos; // Dir from l1 to l2 intersection points
      glm::vec2 dirInv = 1.f/dir;

      // std::println("l1Pos: {}, {}, l2Pos: {}, {}", l1Pos.x, l1Pos.y, l2Pos.x, l2Pos.y);

      drawLine2(dominantAxisIndex, dominantAxisValue, l1Pos, l2Pos, dir, dirInv, pInsertFunc);
    }

    ++*pTrisComplete;
  }
}

void VoxelGrid::writeToFile(const std::string& pPath) {
  std::ofstream fout;
  openFileWrite(fout, pPath);
  writeMetaData(fout);

  fout.write(&mVoxelData.at(0), mVoxelData.size());

  fout.close();
}

void VoxelGrid::writeToFileCompressed(const std::string& pPath, uint* pVoxelsComplete) {
  std::ofstream fout;
  openFileWrite(fout, pPath);
  writeMetaData(fout);

  std::vector<uint> compressedData = generateCompressedVoxelData(pVoxelsComplete);
  fout.write(reinterpret_cast<char*>(&compressedData.at(0)), compressedData.size() * sizeof(uint));

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
  char byte;
  uint i = 0;
  glm::tvec3<uint> pos;
  while (fin.read(&byte, 1)) {
    for (uint bitIndex = 0; bitIndex < 8; ++bitIndex, ++i, byte >>= 1) {
      pos.z = i / (mResolution * mResolution);
      pos.y = (i / mResolution) % mResolution;
      pos.x = i % mResolution;
      if (byte & 1)
        insert(pos);
    }
  }

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

int VoxelGrid::insert(const glm::tvec3<uint>& pPos, const insertFunc_t& pInsertFunc) {
  if (pPos.x < 0 || pPos.x >= mResolution ||
      pPos.y < 0 || pPos.y >= mResolution ||
      pPos.z < 0 || pPos.z >= mResolution) {
    std::lock_guard<std::mutex> lock(*mLogMutex);
    std::println(*mLogStream, "Couldn't insert voxel: {}, {}, {}. Reason: voxel outside grid", pPos.x, pPos.y, pPos.z);
    return 1;
  }

  if (!mVoxelGrid[pPos.x][pPos.y][pPos.z]) {
    ++mVoxelCount;
    pInsertFunc(pPos.x, pPos.y, pPos.z);
    mVoxelGrid[pPos.x][pPos.y][pPos.z] = true;

    uint globalIndex = pPos.x + pPos.y * mResolution + pPos.z * mResolution * mResolution;
    uint byteIndex = globalIndex >> 3;
    uint localIndex = globalIndex & 0b111;
    char mask = 1 << localIndex;
    mVoxelData[byteIndex] |= mask;
  }
  return 0;
}

uint VoxelGrid::getVoxelCount() {
  return mVoxelCount;
}

uint VoxelGrid::getVolume() {
  return mVolume;
}

uint VoxelGrid::getResolution() {
  return mResolution;
}

float VoxelGrid::getMaxDepth() {
  return mMaxDepth;
}

const std::vector<std::vector<std::vector<bool>>>& VoxelGrid::getVoxelData() {
  return mVoxelGrid;
}

const std::vector<char>& VoxelGrid::getVoxelDataBits() {
  return mVoxelData;
}

std::vector<uint> VoxelGrid::generateCompressedVoxelData(uint* pVoxelsComplete) {
  std::vector<uint> counts;
  uint count = 0;
  bool value = false;
  for (uint x = 0; x < mResolution; ++x) {
    for (uint y = 0; y < mResolution; ++y) {
      for (uint z = 0; z < mResolution; ++z) {
        if (mVoxelGrid[x][y][z] != value) {
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
  mVolume = mResolution * mResolution * mResolution;
  mMaxDepth = std::log2f(mResolution);

  mVoxelGrid.clear();
  mVoxelData.clear();

  // Init mVoxelGrid and mVoxelData
  mVoxelGrid = std::vector<std::vector<std::vector<bool>>>(mResolution, std::vector<std::vector<bool>>(mResolution, std::vector<bool>(mResolution, false)));

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

void VoxelGrid::drawLine2(uint8_t pDominantAxisIndex, float pDominantAxisValue, const glm::vec2& pStart, const glm::vec2& pEnd, const glm::vec2& pDir, const glm::vec2& pDirInv, const insertFunc_t& pInsertFunction) {
  glm::vec3 v = glm::floor(toVec3(pDominantAxisValue, pStart.x, pStart.y, pDominantAxisIndex));

  insert(v, pInsertFunction);

  glm::vec2 voxelPos = glm::floor(pStart);

  while (voxelPos != glm::floor(pEnd)) {
    // std::println("voxelPos: {}, {}, l2pos: {}, {}, signDirInv: {}, {}", voxelPos.x, voxelPos.y, pEnd.x, pEnd.y, glm::sign(pDirInv.x), glm::sign(pDirInv.y));
    float plane1 = voxelPos.x + std::max(0.f, glm::sign(pDirInv.x));
    float plane2 = voxelPos.y + std::max(0.f, glm::sign(pDirInv.y));

    if (plane1 == pEnd.x && plane2 == pEnd.y) { // If the destination point lies exactly on the planes we increment both x and y
      voxelPos += glm::sign(pDirInv);
      v = glm::floor(toVec3(pDominantAxisValue, voxelPos.x, voxelPos.y, pDominantAxisIndex));
      insert(v, pInsertFunction);
      break;
    }
    
    float t1 = (plane1 - pStart.x) * pDirInv.x;
    float t2 = (plane2 - pStart.y) * pDirInv.y;

    if (t1 < t2)
      voxelPos.x += glm::sign(pDirInv.x);
    else
      voxelPos.y += glm::sign(pDirInv.y);

    v = glm::floor(toVec3(pDominantAxisValue, voxelPos.x, voxelPos.y, pDominantAxisIndex));

    insert(v, pInsertFunction);
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

