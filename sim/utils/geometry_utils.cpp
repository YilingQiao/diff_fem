#include <iostream>
#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#include <tiny_obj_loader.h>
#undef TINYOBJLOADER_IMPLEMENTATION

#include <utils/geometry_utils.h>


#define DIMENSION 3
#define NB_VERTICES_IN_TRIANGLE 3

bool
read_obj(const std::string& filename,
         std::vector<Eigen::Vector3d>& positions,
         std::vector<tinyobj::index_t>& triangles,
         std::vector<Eigen::Vector3d>& normals,
         std::vector<Eigen::Vector2d>& texcoords)
{
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string err;

    bool ret =
      tinyobj::LoadObj(&attrib, &shapes, &materials, &err, filename.c_str());

    if (!err.empty()) {
        std::cerr << err << std::endl;
    } 

    if (!ret) {
        return ret;
    }

    positions.clear();
    triangles.clear();
    normals.clear();
    texcoords.clear();

    for (int i = 0; i < shapes.size(); i++) {
        assert(shapes[i].mesh.indices.size() % NB_VERTICES_IN_TRIANGLE == 0);
        for (int f = 0; f < shapes[i].mesh.indices.size(); f++) {
            triangles.push_back(shapes[i].mesh.indices[f]);
        }
        assert(attrib.vertices.size() % DIMENSION == 0);
        for (int v = 0; v < attrib.vertices.size() / DIMENSION; v++) {
            positions.push_back(
              Eigen::Vector3d(attrib.vertices[DIMENSION * v + 0],
                              attrib.vertices[DIMENSION * v + 1],
                              attrib.vertices[DIMENSION * v + 2]));
        }

        assert(attrib.normals.size() % DIMENSION == 0);
        for (int n = 0; n < attrib.normals.size() / 3; n++) {
            normals.push_back(Eigen::Vector3d(attrib.normals[3 * n + 0],
                                              attrib.normals[3 * n + 1],
                                              attrib.normals[3 * n + 2]));
        }

        assert(attrib.texcoords.size() % (DIMENSION - 1) == 0);
        for (int t = 0; t < attrib.texcoords.size() / 2; t++) {
            texcoords.push_back(Eigen::Vector2d(attrib.texcoords[2 * t + 0],
                                                attrib.texcoords[2 * t + 1]));
        }
    }

    return ret;
}


#undef DIMENSION
#undef NB_VERTICES_IN_TRIANGLE
