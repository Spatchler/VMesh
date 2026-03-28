#include "palette.hpp"

using namespace VMesh;

Palette::Palette() {
}

uint16_t Palette::addColour(const glm::vec3& pCol, float pThreshold) {
  if (pCol.x > 1 || pCol.x < 0 ||
      pCol.y > 1 || pCol.y < 0 ||
      pCol.z > 1 || pCol.z < 0)
    throw std::runtime_error("Cannot add colour, rgb values have to be between 0 and 1 inclusive");
  for (uint i = 0; i < mColours.size(); ++i) {
    glm::vec3 diff = glm::abs(pCol - mColours[i]);
    if ((diff.x + diff.y + diff.z) / 3.f <= pThreshold) return i;
  }
  mColours.push_back(pCol);
  return mColours.size() - 1;
}

uint16_t Palette::getClosestColour(const glm::vec3& pCol) {
  float minDiff = 1.f;
  uint16_t minDiffIndex = 0;
  for (uint i = 0; i < mColours.size(); ++i) {
    glm::vec3 diff = glm::abs(pCol - mColours[i]);
    float diffPercent = (diff.x + diff.y + diff.z) / 3.f;
    if (diffPercent < minDiff) {
      minDiff = diffPercent;
      minDiffIndex = i;
    }
  }
  return minDiffIndex;
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
  if (!fin.is_open()) throw std::runtime_error("Could not open input palette file");
  // Read
  uint numColours;
  std::string line;
  std::getline(fin, line);
  if (line != "JASC-PAL") throw std::runtime_error("Palette file incorrect format");
  std::getline(fin, line);
  if (line != "0100") throw std::runtime_error("Palette file incorrect format version");
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
      v3index(c, i) = std::stoi(componentStr) / 255.f;
    }
  }
  if (count != numColours) throw std::overflow_error("Palette file invalid number of colours");
  // Close
  fin.close();
}
