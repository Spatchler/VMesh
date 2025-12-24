#include "model.hpp"

using namespace VMesh;

Model::Model(const std::string& pPath) {
  load(pPath);
}

void Model::load(const std::string& pPath) {
  Assimp::Importer importer;

  const aiScene* scene = importer.ReadFile(pPath, aiProcess_Triangulate);

  if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
    throw std::format("ERROR::ASSIMP: {}", importer.GetErrorString());

  processNode(scene->mRootNode, scene);

  mTriCount = mIndices.size() / 3u;
}

void Model::processNode(aiNode* pNode, const aiScene* pScene) {
  for (uint i = 0; i < pNode->mNumMeshes; ++i) {
    aiMesh* mesh = pScene->mMeshes[pNode->mMeshes[i]];
    processMesh(mesh, pScene);
  }

  for (uint i = 0; i < pNode->mNumChildren; ++i)
    processNode(pNode->mChildren[i], pScene);
}

void Model::processMesh(aiMesh* pMesh, const aiScene* pScene) {
  uint numVertsBefore = mVerts.size();

  for (uint i = 0; i < pMesh->mNumVertices; ++i) {
    glm::vec3 pos;

    pos.x = pMesh->mVertices[i].x;
    pos.y = pMesh->mVertices[i].y;
    pos.z = pMesh->mVertices[i].z;

    mVerts.push_back(pos);
  }

  for (uint i = 0; i < pMesh->mNumFaces; ++i) {
    aiFace face = pMesh->mFaces[i];
    for (uint j = 0; j < face.mNumIndices; ++j)
      mIndices.push_back(face.mIndices[j] + numVertsBefore);
  }
}
