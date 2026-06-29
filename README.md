# VMesh

A C++ library for CPU DDA based triangle mesh voxelization. I also plan to make a GPU version with OpenCL.

### Dependencies:
* glm
* assimp
* stb/stb_image
* ogt_vox.h from [opengametools](https://github.com/jpaver/opengametools)

Install dependencies with your system package manager e.g. `# pacman -S glm stb assimp` and they should work straight away since the premake config includes from `/usr/include` and shared objects should be linked by the os

Place `ogt_vox.h` in `dependencies/include`

### Build:
`premake5 gmake && make config=dist`

