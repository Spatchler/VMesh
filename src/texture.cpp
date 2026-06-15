#include "texture.hpp"

using namespace VMesh;

Texture::Texture(const std::string& pPath)
:mPath(pPath), mIsColour(false) {
  std::replace(mPath.begin(), mPath.end(), '\\', '/');
  mData = stbi_load(mPath.c_str(), &mSize.x, &mSize.y, &mNrComponents, 4); // Force 4 channels
  mArea = mSize.x * mSize.y;
  if (!mData) {
    std::println("Texture failed to load at path: {0}", mPath);
    stbi_image_free(mData);
  }
}

Texture::Texture(const glm::vec4& pColour)
:mIsColour(true), mColour(pColour) {
}

glm::vec4 Texture::sample(const glm::vec2& pTexCoord) {
  if (mIsColour) return mColour;
  // std::println("pTexCoord: {}, {}", pTexCoord.x, pTexCoord.y);
  glm::ivec2 pixelCoord = pTexCoord * (glm::vec2(mSize) - glm::vec2(1));
  uint i = ((pixelCoord.y * mSize.x + pixelCoord.x) % mArea) * 4; // 4 because we forced 4 channels not mNrComponents
  uint8_t r = mData[i], g = mData[++i], b = mData[++i], a = mData[++i];
  return glm::vec4(r, g, b, a) / 255.f;
}

const std::string& Texture::getPath() {
  return mPath;
}

void Texture::release() {
  if (!mIsColour) stbi_image_free(mData);
}
