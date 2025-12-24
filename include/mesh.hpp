#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
// #include <glm/gtc/type_ptr.hpp>

#include "timer.hpp"

#include <string>

namespace VMesh {
  class Mesh {
  public:
    Mesh() = default;

    void transformVertices(const glm::mat4& pM);

    void addVert(const glm::vec3& pVert);
    void addTri(const std::array<uint, 3>& pIndices);
    void addTriUnique(const std::array<glm::vec3, 3>& pVerts);

    uint getTriCount() const;

    const std::vector<glm::vec3>& getVertices();
    const std::vector<uint>& getIndices();
  protected:
    std::vector<glm::vec3> mVerts;
    std::vector<uint> mIndices;

    uint mTriCount = 0;
  };
}
