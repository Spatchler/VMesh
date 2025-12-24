#pragma once

#include "mesh.hpp"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

namespace VMesh {
  class Model: public Mesh {
  public:
    Model() = default;
    Model(const std::string& pPath);
    void load(const std::string& pPath);
  private:
    void processNode(aiNode* pNode, const aiScene* pScene);
    void processMesh(aiMesh* pMesh, const aiScene* pScene);
  };
}
