#include "palette.hpp"

using namespace VMesh;

Palette::Palette() {
}

uint8_t Palette::addColour(const glm::vec3& pCol, float pDistance2) {
  if (pCol.x > 1.f || pCol.x < 0.f ||
      pCol.y > 1.f || pCol.y < 0.f ||
      pCol.z > 1.f || pCol.z < 0.f)
    throw std::invalid_argument("rgb values have to be normalized");

  if (pDistance2 != -1.f) for (uint i = 0; i < mColours.size(); ++i)
    if (glm::distance2(pCol, mColours[i]) <= pDistance2) return i;

  mColours.push_back(pCol);
  return mColours.size() - 1;
}

uint8_t Palette::getClosestColour(const glm::vec3& pCol) {
  float minDist2 = 1.f;
  uint8_t minDist2Index = 0;
  for (uint i = 0; i < mColours.size(); ++i) {
    float dist2 = glm::distance2(pCol, mColours[i]);
    if (dist2 < minDist2) {
      minDist2 = dist2;
      minDist2Index = i;
    }
  }
  return minDist2Index;
}

glm::vec3 Palette::getColour(uint8_t pIndex) {
  return mColours.at(pIndex);
}

uint Palette::size() {
  return mColours.size();
}

void Palette::writeToFile(const std::string& pPath) {
  // Open
  std::ofstream fout;
  fout.open(pPath, std::ios::out | std::ios::binary);
  if (!fout.is_open()) throw std::runtime_error("Could not open palette output file");
  // Write
  uint numColours = mColours.size();
  fout << "JASC-PAL\n0100\n" << numColours << "\n";
  for (glm::vec3& c: mColours)
    fout << static_cast<uint>(c.x * 255) << " " << static_cast<uint>(c.y * 255) << " " << static_cast<uint>(c.z * 255) << "\n";
  // Close
  fout.close();
}

void Palette::readFromFile(const std::string& pPath) {
  // Open
  std::ifstream fin;
  fin.open(pPath, std::ios::binary | std::ios::in);
  if (!fin.is_open()) throw std::invalid_argument("Could not open input palette file");
  // Read
  uint numColours;
  std::string line;
  std::getline(fin, line);
  if (line != "JASC-PAL") throw std::invalid_argument("Palette file incorrect format");
  std::getline(fin, line);
  if (line != "0100") throw std::invalid_argument("Palette file incorrect format version");
  std::getline(fin, line);
  numColours = std::stoi(line);
  mColours.resize(numColours);
  uint count = 0;
  for (; std::getline(fin, line) && (count < numColours); ++count) {
    std::stringstream s;
    s << line;
    glm::vec3& c = mColours.at(count);
    std::string componentStr;
    for (uint i = 0; i < 3; ++i) {
      std::getline(s, componentStr, ' ');
      c[i] = std::stoi(componentStr) / 255.f;
    }
  }
  if (count != numColours) throw std::overflow_error("Palette file invalid number of colours");
  // Close
  fin.close();
}
