#pragma once

#include <stb_image.h>

#include <glm/glm.hpp>

#include <string>
#include <algorithm>
#include <print>

namespace VMesh {
  class Texture {
  public:
    Texture(const std::string& pPath);
    Texture(const glm::vec4& pColour);

    glm::vec4 sample(const glm::vec2& pTexCoord);

    const std::string& getPath();

    void release();
  private:
    bool mIsColour;
    glm::vec4 mColour;
    std::string mPath;
    uint8_t* mData;
    glm::ivec2 mSize;
    uint mArea;
    int mNrComponents;
  };
}
