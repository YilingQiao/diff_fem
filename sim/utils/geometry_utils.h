
#ifndef GEOMETRY_UTILS_H
#define GEOMETRY_UTILS_H
#include <Eigen/Core>
#include <numeric>
#include <string>
#include <tiny_obj_loader.h>
#include <vector>

template <typename TinyScalar, typename TinyConstants> 
Eigen::VectorBlock<Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>, 3>
getVector3dBlock(Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& vector, int index)
{
    return vector.template segment<3>(3 * index);
}


// {
//     return vector.segment<3>(3 * index);
// }
template <typename TinyScalar, typename TinyConstants>
const Eigen::VectorBlock<const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>, 3>
getVector3dBlock(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& vector, int index)
{
    return vector.template segment<3>(3 * index);
}
// {
//     return vector.segment<3>(3 * index);
// }

template <typename TinyScalar, typename TinyConstants>
std::size_t
getVector3dSize(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& vector)
{
    return vector.size() / 3;
}


bool
read_obj(const std::string& filename,
         std::vector<Eigen::Vector3d>& positions,
         std::vector<tinyobj::index_t>& triangles,
         std::vector<Eigen::Vector3d>& normals,
         std::vector<Eigen::Vector2d>& texcoords);

#endif