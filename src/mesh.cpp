#include "mesh.hpp"

using namespace VMesh;

void Mesh::transformVertices(const glm::mat4& pM) {
  for (uint i = 0; i < mVerts.size(); ++i)
    mVerts[i] = glm::vec3(pM * glm::vec4(mVerts[i], 1));
}

void Mesh::addVert(const glm::vec3& pVert) {
  mVerts.push_back(pVert);
}

void Mesh::addTri(const std::array<uint, 3>& pIndices) {
  mIndices.insert(mIndices.end(), pIndices.begin(), pIndices.end());
  ++mTriCount;
}

void Mesh::addTriUnique(const std::array<glm::vec3, 3>& pVerts) {
  mVerts.insert(mVerts.end(), pVerts.begin(), pVerts.end());
  for (int i = 0; i >= -2; --i)
    mIndices.push_back(mVerts.size() - i);
  ++mTriCount;
}

uint Mesh::getTriCount() const {
  return mTriCount;
}

const std::vector<glm::vec3>& Mesh::getVertices() {
  return mVerts;
}

const std::vector<uint>& Mesh::getIndices() {
  return mIndices;
}
