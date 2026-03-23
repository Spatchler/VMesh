#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
// #include <glm/gtc/type_ptr.hpp>

#include "timer.hpp"

#include "texture.hpp"

#include <string>

namespace VMesh {
  struct Vertex {
    glm::vec3 pos;
    glm::vec2 texCoord;
  };

  class Mesh {
  public:
    Mesh() = default;

    void setMatIndex(uint i);
    uint getMatIndex();

    void transformVertices(const glm::mat4& pM);

    void addVertex(const Vertex& pVert);
    void addTri(const std::array<uint, 3>& pIndices);
    void addTriUnique(const std::array<Vertex, 3>& pVerts);

    uint getTriCount() const;

    const std::vector<Vertex>& getVertices();
    const std::vector<uint>& getIndices();
  protected:
    uint mMatIndex;
    std::vector<Vertex> mVerts;
    std::vector<uint> mIndices;

    uint mTriCount = 0;
  };
}
