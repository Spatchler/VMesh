#include "texture.hpp"

using namespace VMesh;

Texture::Texture(const std::string& pPath)
:mPath(pPath), mIsColour(false) {
  std::replace(mPath.begin(), mPath.end(), '\\', '/');
  mData = stbi_load(mPath.c_str(), &mSize.x, &mSize.y, &mNrComponents, 3); // Force 3 channels
  mArea = mSize.x * mSize.y;
  if (!mData) {
    std::println("Texture failed to load at path: {0}", mPath);
    stbi_image_free(mData);
  }
}

Texture::Texture(const glm::vec3& pColour)
:mIsColour(true), mColour(pColour) {
}

glm::vec3 Texture::sample(const glm::vec2& pTexCoord) {
  if (mIsColour) return mColour;
  // std::println("sampling");
  glm::ivec2 pixelCoord = glm::ivec2(glm::abs(pTexCoord)) * mSize;
  uint i = ((pixelCoord.y * mSize.y) + pixelCoord.x) * mNrComponents % mArea;
  uint8_t r = *(mData+i), g = *(mData+i+1), b = *(mData+i+2);
  // std::println();
  // uint8_t r = mData[i], g = mData[++i], b = mData[++i];
  return glm::vec3(r, g, b) / 255.f;
}

const std::string& Texture::getPath() {
  return mPath;
}

void Texture::release() {
  if (!mIsColour) stbi_image_free(mData);
}
