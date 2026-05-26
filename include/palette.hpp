#pragma once

#include <glm/glm.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/norm.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <print>

namespace VMesh {
  class Palette {
  public:
    Palette();

    uint8_t addColour(const glm::vec3& pCol, float pDistance2 = 0.f);
    uint8_t getClosestColour(const glm::vec3& pCol);
    glm::vec3 getColour(uint8_t pIndex);

    uint size();

    void writeToFile(const std::string& pPath);
    void readFromFile(const std::string& pPath);
  private:
    std::vector<glm::vec3> mColours;
  };
}
