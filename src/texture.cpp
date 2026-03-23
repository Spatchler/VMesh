using namespace VMesh;

Texture::Texture(const std::string& pPath)
:mData = stbi_load(pPath.c_str(), &mSize.x, &mSize.y, &mNrComponents, 3), mPath(pPath) {
  if (!mData) {
    std::println("Texture failed to load at path: {0}", pPath);
    stbi_image_free(mData);
  }
}

glm::vec3 Texture::sample(const glm::vec2& pTexCoord) {
  glm::ivec2 pixelCoord = pTexCoord * mSize;
  uint i = (pixelCoord.y / mSize.y * mSize.y + pixelCoord.x) * 3;
  uint8_t r = mData[i], g = mData[++i], b = mData[++i];
  return glm::vec3(r, g, b) / 255;
}

const std::string& Texture::getPath() {
  return mPath;
}

Texture::~Texture() {
  stbi_image_free(mData);
}
