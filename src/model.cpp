#include "model.hpp"

using namespace VMesh;

Model::Model(const std::string& pPath) {
  load(pPath);
}

void Model::load(const std::string& pPath, bool pShouldLoadTexture) {
  Assimp::Importer importer;

  const aiScene* scene = importer.ReadFile(pPath, aiProcess_Triangulate);

  if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
    throw std::runtime_error(importer.GetErrorString());

  mDirectory = pPath.substr(0, pPath.find_last_of('/'));

  loadMaterials(scene);

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

void Model::release() {
  for (auto& i: mDiffuseMaps) i.release();
}

void Model::processNode(aiNode* pNode, const aiScene* pScene) {
  for (uint i = 0; i < pNode->mNumMeshes; ++i)
    processMesh(pScene->mMeshes[pNode->mMeshes[i]], pScene);

  for (uint i = 0; i < pNode->mNumChildren; ++i)
    processNode(pNode->mChildren[i], pScene);
}

void Model::processMesh(aiMesh* pMesh, const aiScene* pScene) {
  mMeshes.emplace_back();
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
    mMeshes.back().addTri({face.mIndices[0], face.mIndices[1], face.mIndices[2]});
    ++mTriCount;
  }
}

void Model::loadMaterials(const aiScene* pScene) {
  if (!pScene->HasMaterials()) return;
  for (uint i = 0; i < pScene->mNumMaterials; ++i) {
    aiMaterial* mat = pScene->mMaterials[i];
    // Base Colour if has no textures
    if (mat->GetTextureCount(aiTextureType_DIFFUSE) == 0) {
      aiColor3D colour(0.f,0.f,0.f);
      mat->Get(AI_MATKEY_COLOR_DIFFUSE, colour);
      mDiffuseMaps.emplace_back(glm::vec3(colour.r, colour.g, colour.b));
      continue;
    }
    // Texture
    for (uint j = 0; j < mat->GetTextureCount(aiTextureType_DIFFUSE); ++j) {
      aiString path;
      mat->GetTexture(aiTextureType_DIFFUSE, j, &path);
      mDiffuseMaps.emplace_back(mDirectory + '/' + std::string(path.C_Str()));
    }
  }
}
