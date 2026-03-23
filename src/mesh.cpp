#include "mesh.hpp"

using namespace VMesh;

void Mesh::setMatIndex(uint i) {
  mMatIndex = i;
}

uint Mesh::getMatIndex() {
  return mMatIndex;
}

void Mesh::transformVertices(const glm::mat4& pM) {
  for (uint i = 0; i < mVerts.size(); ++i)
    mVerts[i] = glm::vec3(pM * glm::vec4(mVerts[i], 1));
}

void Mesh::addVertex(const Vertex& pVert) {
  mVerts.push_back(pVert);
}

void Mesh::addTri(const std::array<uint, 3>& pIndices) {
  mIndices.insert(mIndices.end(), pIndices.begin(), pIndices.end());
  ++mTriCount;
}

void Mesh::addTriUnique(const std::array<Vertex, 3>& pVerts) {
  mVerts.insert(mVerts.end(), pVerts.begin(), pVerts.end());
  for (int i = 0; i >= -2; --i)
    mIndices.push_back(mVerts.size() - i);
  ++mTriCount;
}

uint Mesh::getTriCount() const {
  return mTriCount;
}

const std::vector<Vertex>& Mesh::getVertices() {
  return mVerts;
}

const std::vector<uint>& Mesh::getIndices() {
  return mIndices;
}
