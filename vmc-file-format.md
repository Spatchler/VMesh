### Description:

File format for storing 8-bit voxel data compressed with run length encoding. Often paired with a JASC-PAL file.

### Contens:

| Bytes    | Type       | Value                                                  |
| :------- | :--------- | :----------------------------------------------------- |
| 1\*6     | char       | id 'VMESHC' : 'V' 'M' 'E' 'S' 'H' 'C', 'V' is first    |
| 4        | uint       | version number : 100                                   |
| 4        | uint       | grid resolution                                        |
| 1        | uint       | palette size not including air                         |
| 8        | uint       | voxel count                                            |
| 4        | uint       | (N)                                                    |
| 1\*N     | uint       | Voxel types array                                      |
| 8\*N     | uint       | Voxel counts array                                     |

### Loading and decompression example:

Resulting voxel data is ordered using a z order curve.

```cpp
void load(const std::filesystem::path& path, uint32_t& resolution, uint8_t& paletteSize, uint64_t& voxelCount, std::vector<uint8_t>& voxelData) {
  // Open file
  std::ifstream fin;
  fin.open(path, std::ios::binary | std::ios::in);
  if (!fin.is_open()) throw std::runtime_error("Could not open file");
  
  // Check header
  std::string str;
  str.resize(6);
  fin.read(str.data(), 6);
  if (str != "VMESHC") throw std::invalid_argument("Invalid voxel grid file");
  uint32_t version;
  fin.read(reinterpret_cast<char*>(&version), sizeof(version));
  if (ver != 100) throw std::invalid_argument("Invalid voxel grid file version");
  
  // Read metadata
  fin.read(reinterpret_cast<char*>(&resolution), sizeof(uint32_t));
  fin.read(reinterpret_cast<char*>(&paletteSize), sizeof(uint8_t));
  fin.read(reinterpret_cast<char*>(&voxelCount), sizeof(uint64_t));
  
  // Read counts len
  uint32_t countsLen;
  fin.read(reinterpret_cast<char*>(&countsLen), sizeof(uint32_t));
  voxelData.resize(volume, 0);
  
  // Read voxel data
  std::vector<uint8_t> values;
  values.resize(countsLen);
  fin.read(reinterpret_cast<char*>(values.data()), countsLen * sizeof(uint8_t));
  
  std::vector<uint64_t> counts;
  counts.resize(countsLen);
  fin.read(reinterpret_cast<char*>(counts.data()), countsLen * sizeof(uint64_t));
  
  // Decompress
  uint64_t volume = static_cast<uint64_t>(resolution) * resolution * resolution;
  voxelData.resize(volume, 0);
  uint64_t index = 0;
  for (uint i = 0; i < countsLen; ++i) {
    const uint8_t& value = values[i];
    const uint64_t& count = counts[i];
    if (value != 0) for (uint64_t j = 0; j < count; ++j) voxelData[index + j] = value;
    index += count;
  }
  
  // Close
  fin.close();
}
```

### Z-Order curve function example:

source: [Bit Twiddling Hacks](https://graphics.stanford.edu/%7Eseander/bithacks.html#InterleaveBMN)
```cpp
uint64_t interleave(uint64_t x, uint64_t y, uint64_t z) {
  static const uint64_t B[] = {0x00000000FF0000FF, 0x000000F00F00F00F,
                               0x00000C30C30C30C3, 0X0000249249249249};
  static const int S[] =  {16, 8, 4, 2}; 

  x = (x | (x << S[0])) & B[0];
  x = (x | (x << S[1])) & B[1];
  x = (x | (x << S[2])) & B[2];
  x = (x | (x << S[3])) & B[3];

  y = (y | (y << S[0])) & B[0];
  y = (y | (y << S[1])) & B[1];
  y = (y | (y << S[2])) & B[2];
  y = (y | (y << S[3])) & B[3];

  z = (z | (z <<  S[0])) & B[0];
  z = (z | (z <<  S[1])) & B[1];
  z = (z | (z <<  S[2])) & B[2];
  z = (z | (z <<  S[3])) & B[3];

  return x | (y << 1) | (z << 2);
}
```
