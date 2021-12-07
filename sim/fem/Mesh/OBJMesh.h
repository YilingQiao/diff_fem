#ifndef __OBJ_MESH_H__
#define __OBJ_MESH_H__
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <fstream>
#include <sstream>
#include <iostream>
#include "Mesh.h"

namespace FEM
{
template <typename TinyScalar, typename TinyConstants> 
class OBJMesh : public Mesh<TinyScalar, TinyConstants>
{
public:
	OBJMesh(const std::string& obj_file,const Eigen::Transform<TinyScalar, 3, 2>& T = Eigen::Transform<TinyScalar, 3, 2>::Identity());	
};




using namespace FEM;

template <typename TinyScalar, typename TinyConstants> 
OBJMesh<TinyScalar, TinyConstants>::
OBJMesh(const std::string& path,const Eigen::Transform<TinyScalar, 3, 2>& T)
    :Mesh<TinyScalar, TinyConstants>()
{
    std::ifstream ifs(path);
    if(!(ifs.is_open()))
    {
        std::cout<<"Can't read file "<<path<<std::endl;
        return;
    }
    std::string str;
    std::string index;
    std::stringstream ss;

    while(!ifs.eof())
    {
        str.clear();
        index.clear();
        ss.clear();

        std::getline(ifs,str);
        ss.str(str);
        ss>>index;

        if(!index.compare("v"))
        {
            TinyScalar x,y,z;
            ss>>x>>y>>z;
            this->mVertices.push_back(Eigen::Matrix<TinyScalar, 3, 1>(x,y,z));
        }
        else if(!index.compare("vt"))
        {
            TinyScalar x,y;
            ss>>x>>y;
            this->mTextureCoord.push_back(Eigen::Matrix<TinyScalar, 2, 1>(x,y));
        }
        else if(!index.compare("vn"))
        {
            TinyScalar x,y,z;
            ss>>x>>y>>z;
            this->mVerticesNormal.push_back(Eigen::Matrix<TinyScalar, 3, 1>(x,y,z));
        }
        else if(!index.compare("f"))
        {
            int i0,i1,i2;
            int n0,n1,n2;
            int t0, t1, t2;
            Eigen::Matrix<TinyScalar, 3, 1> position;
            Eigen::Matrix<TinyScalar, 3, 1> normal;

            const char* chh=str.c_str();
            sscanf (chh, "f %d/%d/%d %d/%d/%d %d/%d/%dm",
                &i0,&t0,&n0,
                &i1,&t1,&n1,
                &i2,&t2,&n2);

            this->mTriangles.push_back(Eigen::Vector3i(i0,i1,i2));
            this->mFacesNormal.push_back(Eigen::Vector3i(n0,n1,n2));
            this->mFacesTexture.push_back(Eigen::Vector3i(t0,t1,t2));
        }
        else if(!index.compare("t"))
        {
            int i0,i1,i2,i3;
            ss>>i0>>i1>>i2>>i3;
            this->mTetrahedrons.push_back(Eigen::Vector4i(i0,i1,i2,i3));
        }

    }
    ifs.close();

    for(auto& v : this->mVertices)
        v = T*v;
}

};
#endif