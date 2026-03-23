#pragma once

#include <stb_image.h>

#include <string>

namespace VMesh {
  class Texture {
  public:
    Texture(const std::string& pPath);

    glm::vec3 sample(const glm::vec2& pTexCoord);

    const std::string& getPath();

    ~Texture();
  private:
    const std::string mPath;
    uint8_t* mData;
    glm::vec2 mSize;
    int mNrComponents;
  };
}
