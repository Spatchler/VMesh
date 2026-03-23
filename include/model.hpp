#pragma once

#include "mesh.hpp"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

namespace VMesh {
  class Model {
  public:
    Model() = default;
    Model(const std::string& pPath);
    void load(const std::string& pPath, bool pShouldLoadTexture = false);
    Texture& getDiffuseMap(uint i);
    uint getTriCount();
    uint getNumMeshes();
    Mesh& getMesh(uint i);
  private:
    uint mTriCount = 0;
    std::vector<Mesh> mMeshes;
    std::string mDirectory;
    std::vector<Texture> mDiffuseMaps;
    void loadMaterialDiffuseMaps(aiMaterial* pMat);
    void processNode(aiNode* pNode, const aiScene* pScene);
    void processMesh(aiMesh* pMesh, const aiScene* pScene);
  };
}
