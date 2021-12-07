#ifndef __OBJ_CLOTH_MESH_H__
#define __OBJ_CLOTH_MESH_H__
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <fstream>
#include <sstream>
#include <iostream>
#include "Mesh.h"
#include "CollisionMesh.h"
#include "../../utils/geometry_utils.h"

namespace FEM
{

template <typename TinyScalar, typename TinyConstants> 
class OBJClothMesh : public Mesh<TinyScalar, TinyConstants>
{
public:
    OBJClothMesh(const std::string& obj_file,const Eigen::Transform<TinyScalar, 3, 2>& T = Eigen::Transform<TinyScalar, 3, 2>::Identity());    
    std::vector<tinyobj::index_t> mTinyTriangles;
    ClothPhysicParameters<TinyScalar, TinyConstants> mParameter;
};



static void get_valid_line (std::istream &in, std::string &line) {
    do
        getline(in, line);
    while (in && (line.length() == 0 || line[0] == '#'));
}

using namespace FEM;

template <typename TinyScalar, typename TinyConstants> 
OBJClothMesh<TinyScalar, TinyConstants>::
OBJClothMesh(const std::string& path,const Eigen::Transform<TinyScalar, 3, 2>& T)
    :Mesh<TinyScalar, TinyConstants>()
{

    std::cout << "read obj " << path << std::endl;


    std::vector<Eigen::Vector3d> tmp_positions;
    std::vector<Eigen::Vector3d> tmp_normals;
    std::vector<Eigen::Vector2d> tmp_texcoords;

    bool success;
    success = read_obj(path, tmp_positions, 
    this->mTinyTriangles, tmp_normals, tmp_texcoords);

    this->mVertices = 
        helper::to_eigen_double2tiny<TinyScalar, TinyConstants>(tmp_positions);
    this->mVerticesNormal = 
        helper::to_eigen_double2tiny<TinyScalar, TinyConstants>(tmp_normals);
    this->mTextureCoord = 
        helper::to_eigen_double2tiny<TinyScalar, TinyConstants>(tmp_texcoords);

    for (std::size_t i = 0; i < this->mTinyTriangles.size() / 3; ++i)
    {
        int f0, f1, f2;
        f0 = this->mTinyTriangles[3 * i].vertex_index;
        f1 = this->mTinyTriangles[3 * i + 1].vertex_index;
        f2 = this->mTinyTriangles[3 * i + 2].vertex_index;
        this->mTriangles.push_back(Eigen::Vector3i(f0, f1, f2));
        // TODO: Auxiliary function to get triangle edges
    } // min face number is zero
}

};
#endif