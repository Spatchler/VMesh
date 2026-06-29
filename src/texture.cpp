#include "texture.hpp"

using namespace VMesh;

Texture::Texture(const std::string& pPath, const std::array<aiTextureMapMode, 3>& pMapModes)
:mPath(pPath), mIsColour(false), mMapModes(pMapModes) {
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

glm::vec4 Texture::sample(glm::vec2 pTexCoord) {
  if (mIsColour) return mColour;
  if ((pTexCoord.x > 1.f) || (pTexCoord.x < 0.f)) {
    if (mMapModes[0] == aiTextureMapMode_Wrap) pTexCoord.x = pTexCoord.x - std::floor(pTexCoord.x);
    else if (mMapModes[0] == aiTextureMapMode_Clamp) pTexCoord.x = std::clamp(pTexCoord.x, 0.f, 1.f);
    else if (mMapModes[0] == aiTextureMapMode_Mirror) {
      bool isOddTile = (static_cast<int>(std::floor(pTexCoord.x)) & 1) != 0;
      pTexCoord.x = isOddTile ? (1.f - (pTexCoord.x - std::floor(pTexCoord.x))) : (pTexCoord.x - std::floor(pTexCoord.x));
    }
  }
  if ((pTexCoord.y > 1.f) || (pTexCoord.y < 0.f)) {
    if (mMapModes[1] == aiTextureMapMode_Wrap) pTexCoord.y = pTexCoord.y - std::floor(pTexCoord.y);
    else if (mMapModes[1] == aiTextureMapMode_Clamp) pTexCoord.y = std::clamp(pTexCoord.y, 0.f, 1.f);
    else if (mMapModes[1] == aiTextureMapMode_Mirror) {
      bool isOddTile = (static_cast<int>(std::floor(pTexCoord.y)) & 1) != 0;
      pTexCoord.y = isOddTile ? (1.f - (pTexCoord.y - std::floor(pTexCoord.y))) : (pTexCoord.y - std::floor(pTexCoord.y));
    }
  }
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
