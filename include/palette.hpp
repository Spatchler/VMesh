#pragma once

#include <glm/glm.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <print>

namespace VMesh {
  constexpr float& v3index(glm::tvec3<float>& v, uint8_t i) {
    if (i == 0) return v.x;
    if (i == 1) return v.y;
    if (i == 2) return v.z;
    return v.x;
  }

  class Palette {
  public:
    Palette();

    uint16_t addColour(const glm::vec3& pCol, float pThreshold = 0.f);
    uint16_t getClosestColour(const glm::vec3& pCol);

    void writeToFile(const std::string& pPath);
    void readFromFile(const std::string& pPath);
  private:
    std::vector<glm::vec3> mColours;
  };
}
