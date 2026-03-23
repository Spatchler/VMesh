#include "model.hpp"

using namespace VMesh;

Model::Model(const std::string& pPath) {
  load(pPath);
}

void Model::load(const std::string& pPath, bool pShouldLoadTexture) {
  Assimp::Importer importer;

  const aiScene* scene = importer.ReadFile(pPath, aiProcess_Triangulate);

  if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
    throw std::format("ERROR::ASSIMP: {}", importer.GetErrorString());

  mDirectory = pPath.substr(0, pPath.find_last_of('/'));

  loadMaterials();

  processNode(scene->mRootNode, scene);
}

Texture& Model::getDiffuseMap(uint i) {
  return mDiffuseMaps[i];
}

uint Model::getTriCount() {
  return mTriCount;
}

uint Model::getNumMeshes() {
  return mMeshes.size();
}

Mesh& Model::getMesh(uint i) {
  return mMeshes.at(i);
}

void Model::processNode(aiNode* pNode, const aiScene* pScene) {
  for (uint i = 0; i < pNode->mNumMeshes; ++i)
    processMesh(pScene->mMeshes[pNode->mMeshes[i]], pScene);

  for (uint i = 0; i < pNode->mNumChildren; ++i)
    processNode(pNode->mChildren[i], pScene);
}

void Model::processMesh(aiMesh* pMesh, const aiScene* pScene) {
  mMeshes.empace_back();
  mMeshes.back().setMatIndex(pMesh->mMaterialIndex);
  for (uint i = 0; i < pMesh->mNumVertices; ++i) {
    Vertex v;

    v.pos.x = pMesh->mVertices[i].x;
    v.pos.y = pMesh->mVertices[i].y;
    v.pos.z = pMesh->mVertices[i].z;

    if (pMesh->mTextureCoords[0]) {
      v.texCoord.x = pMesh->mTextureCoords[0][i].x;
      v.texCoord.y = pMesh->mTextureCoords[0][i].y;
    }
    else
      v.texCoord = glm::vec2(0.f, 0.f);

    mMeshes.back().addVertex(v);
  }

  for (uint i = 0; i < pMesh->mNumFaces; ++i) {
    aiFace face = pMesh->mFaces[i]; // A face will always have 3 vertices because we use aiProcess_Triangulate
    mMeshes.back().addTri(face.mIndices[0], face.mIndices[1], face.mIndices[2]);
    ++mTriCount;
  }
}

void Model::loadMaterials(aiScene* pScene) {
  if (!pScene->HasMaterials()) return;
  for (aiMaterial* mat: pScene->mMaterials) {
    for (uint i = 0; i < pMat->GetTextureCount(aiTextureType_DIFFUSE); ++i) {
      aiString path;
      pMat->GetTexture(aiTextureType_DIFFUSE, i, &path);
      mDiffuseMaps.emplace_back(mDirectory + '/' + std::string(path.C_Str()));
    }
  }
}
