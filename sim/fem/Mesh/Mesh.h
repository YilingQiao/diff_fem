#ifndef __MESH_H__
#define __MESH_H__
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <memory>
namespace FEM
{
template <typename TinyScalar, typename TinyConstants> 
class Mesh
{
public:
	Mesh(){};
	virtual const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& GetVertices(){return mVertices;};
	virtual const std::vector<Eigen::Vector3i>& GetTriangles(){return mTriangles;};
	virtual const std::vector<Eigen::Vector4i>& GetTetrahedrons(){return mTetrahedrons;};
	virtual const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& GetVertexNormal() {return mVerticesNormal;};
	virtual const std::vector<Eigen::Matrix<TinyScalar, 2, 1>>& GetTextureCoord() {return mTextureCoord;};
	virtual const std::vector<Eigen::Vector3i>& GetFaceNormal() {return mFacesNormal;};
	virtual const std::vector<Eigen::Vector3i>& GetFaceTexture() {return mFacesTexture;};
	virtual const std::vector<Eigen::Vector2i>& GetFaceEdges() {return mEdges;};
	virtual void Clear() {mVertices.clear(); mTetrahedrons.clear();};
	
// protected:
	std::vector<std::vector<int>> mAdjacentVerts;
	std::vector<Eigen::Matrix<TinyScalar, 3, 1>> mVertices;
	std::vector<Eigen::Vector3i> mTriangles;
	std::vector<Eigen::Vector4i> mTetrahedrons;
	std::vector<Eigen::Matrix<TinyScalar, 3, 1>> mVerticesNormal;
	std::vector<Eigen::Matrix<TinyScalar, 2, 1>> mTextureCoord;
	std::vector<Eigen::Vector3i> mFacesNormal;
	std::vector<Eigen::Vector3i> mFacesTexture;
	std::vector<Eigen::Vector2i> mEdges;



	TinyScalar getInitialEdgeLength(const Eigen::Vector2i& edge) const
	{

	    return (mVertices[edge[0]]-mVertices[edge[1]]).norm();
	}

	std::vector<Eigen::Matrix<TinyScalar, 3, 1>> 
	getInitialVertices(const std::vector<int>& vertex_indices) const
	{
		std::vector<Eigen::Matrix<TinyScalar, 3, 1>> results;
		for (const auto& vi : vertex_indices)
			results.push_back(mVertices[vi]);
		return results;
	}
};


};


#endif