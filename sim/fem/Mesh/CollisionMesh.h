/*
Copyright 2021 by Inria, Mickaël Ly, Jean Jouve, Florence Bertails-Descoubes and
    Laurence Boissieux

This file is part of ProjectiveFriction.

ProjectiveFriction is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

ProjectiveFriction is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
ProjectiveFriction. If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef COLLISION_MESH_H
#define COLLISION_MESH_H

#include <CGAL/Surface_mesh.h>
#include <boost/iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/sub_range.hpp>
#include <Eigen/Core>
#include <functional>
#include <tiny_obj_loader.h>
#include <vector>
#include <fstream>

#include "../Constraint/ConstraintHeader.h"
#include "../../utils/eigen_help.h"
#include "../../utils/geometry_utils.h"
#include "../../diff/tiny_double_utils.h"
#include "../../diff/tiny_eigen_helper.h"

/** @class CollisionMesh
 * A surface mesh made of triangles.
 */

namespace FEM
{


template <typename TinyScalar, typename TinyConstants> 
struct ClothPhysicParameters
{
  /// @brief Bending stiffness
  TinyScalar bend;
  /// @brief Stretch stiffness
  TinyScalar stretch;
  /// @brief Mass per surface unit
  TinyScalar area_density;
};


template <typename TinyScalar, typename TinyConstants> 
class CollisionMesh
{
public:
    // typedef CGAL::Surface_mesh<Eigen::Matrix<TinyScalar, 3, 1> > SurfaceMesh;
    using SurfaceMesh = CGAL::Surface_mesh<Eigen::Matrix<TinyScalar, 3, 1> >;
    // using SurfaceMesh = CGAL::Surface_mesh<Eigen::Vector3d >;
  
    /**
     * Type of an iterator over the interior vertices of the mesh.
     * @see getInteriorVerticesEnd
     * @see getInteriorVerticesBegin
     */
    using InteriorVerticesIterator = boost::filter_iterator<
        std::function<bool(std::size_t)>,
        typename SurfaceMesh::Vertex_iterator>;
    /**
     * Type of a range over the indices of vertices adjacent to a vertex.
     * @see getVerticesAroundVertexRange
     */
    using VerticesAroundVertexRange =
      CGAL::Iterator_range<CGAL::Vertex_around_target_iterator<SurfaceMesh>>;
    /**
     * Type of a range over the indices of the mesh's vertices.
     */
    using VerticesRange = typename SurfaceMesh::Vertex_range;
    using Vertex = Eigen::Matrix<TinyScalar, 3, 1>; /**< Type of the vertices position \see getVertex */

    /**
     * An edge of the mesh.
     */
    struct Edge
    {
        /** Stores the indices of the vertices the edge is made of. */
        std::size_t vertex1_index, vertex2_index;
    };

    /**
     * The type of an iterator over every edge within the mesh.
     * @see getEdgeBegin
     * @see getEdgeEnd
     */
    using EdgeIterator = boost::transform_iterator<
      std::function<Edge(const typename SurfaceMesh::Edge_index&)>,
      typename SurfaceMesh::Edge_iterator>;
    /**
     * The type of a range over every edge within the mesh. 
     * @see getEdgeRange
     */
    using EdgeRange = boost::iterator_range<EdgeIterator>;

    /**
     * A triangle of the mesh.
     */
    struct Triangle
    {
        /** Stores the indices of the vertices the triangle is made of */
        std::array<int, 3> vertex_indices;
    };
    /**
     * The type of an iterator over every triangle within the mesh.
     * @see getTriangleBegin
     * @see getTriangleEnd
     */
    using AllTriangleIterator = boost::transform_iterator<
      std::function<Triangle(const typename SurfaceMesh::Face_index&)>,
        typename SurfaceMesh::Face_iterator>;
    /**
     * The type of a range over every triangle within the mesh.
     * @see getTriangleRange
     */
    using AllTriangleRange = boost::iterator_range<AllTriangleIterator>;
    using TriangleIterator = boost::transform_iterator<
      std::function<Triangle(const typename SurfaceMesh::Face_index&)>,
      boost::filter_iterator<
        std::function<bool(const typename SurfaceMesh::Face_index&)>,
        typename SurfaceMesh::Face_around_target_iterator>>;
    using TriangleRange = boost::iterator_range<TriangleIterator>;

    CollisionMesh(const CollisionMesh& other) = default;
    CollisionMesh(CollisionMesh&& other) = default;

    /**
     * Contruct a mesh from a set of ordered vertices and a set of triangles.
     *
     * @param positions Vector of orderd vertices. The order is used to know
     *                  which vertices are in a triangle.
     *
     * @param triangles Vector of triangle description following
     *                  tiny_obj_reader usage.
     *
     * @see tiny_obj_reader.h
     */
    CollisionMesh(const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& positions,
         const std::vector<int>& triangles);
    CollisionMesh(const SurfaceMesh& surface_mesh);

    /**
     * Apply a tranformation on each vertices of the mesh.
     * @param transform The transformation to apply.
     */
    template<int Mode, int Option>
    void applyTransformation(const Eigen::Transform<TinyScalar, 3, Mode, Option>& transform) noexcept;

    /**
     * Returns the CGAL Surface Mesh associated to the mesh. This method is
     * mostly used for outputing the data in a file. The returned mesh might
     * not be storing the up to date position of the vertices, you might want to
     * call the method updateUnderlyingMesh before calling this method.
     */
    const SurfaceMesh& getUnderlyingMesh() const;
    /**
     * Updates the position of the vertices stored in the CGAL Mesh. This method
     * is mostly used for outputing data in a file in conjunction to
     * getUnderlyingMesh.
     */
    void updateUnderlyingMesh();

    /**
     * Returns the number of vertices in the mesh.
     */
    std::size_t getNumberOfVertices() const;
    /**
     * Returns the number of triangles in the mesh.
     */
    std::size_t getNumberOfTriangles() const noexcept;

    /**
     * Returns the area of the given triangle.
     */
    TinyScalar getTriangleArea(const Triangle& triangle) const;
    /**
     * Returns the outward normal of the given triangle. The normal is obtained
     * through the normalized cross product of the edge (v0, v1) and the edge
     * (v0, v2) of the triangle (in this order). 
     */
    Eigen::Matrix<TinyScalar, 3, 1> getTriangleNormal(const Triangle& triangle) const;
    /**
     * Returns the axis aligned bounding box of the given triangle.
     */
    CGAL::Bbox_3 getTriangleBoundingBox(const Triangle& triangle) const;
    /**
     *  @brief Return the position of the triangle vertices in the column of a
     *  matrix.  Denoting the result by `result`, `result.col(i)` is the
     *  position of the `i`th vertex of the triangle.
     */
    Eigen::Matrix<TinyScalar, 3, 3> getTriangle(const Triangle& triangle) const noexcept;
    /**
     * Returns the triangle which has the given index. 
     */
    CollisionMesh::Triangle getTriangle(std::size_t triangle_index) const noexcept;

    /**
     * Returns a range over every triangle of the mesh. The returned range
     * iterates over object of type Triangle.
     * @getTriangleBegin
     * @getTriangleEnd
     */
    AllTriangleRange getTriangleRange() const noexcept;
    /**
     * Returns the starting iterator of a range over every triangles in the
     * mesh. This iterator iterates on object of type Triangle.
     * @see getTriangleRange
     * @see getTriangleEnd
     */
    AllTriangleIterator getTriangleBegin() const noexcept;
    /**
     * Returns the after end iterator of a range over every triangles in the
     * mesh.
     * @see getTriangleRange
     * @see getTriangleEnd
     */
    AllTriangleIterator getTriangleEnd() const noexcept;

    /**
     * Returns a range over every triangle incident to the vertex that has the
     * given index within the mesh. The returned range iterates over object of
     * type Triangle.
     * @see getTriangleAroundVertexBegin
     * @see getTriangleAroundVertexEnd
     */
    TriangleRange getTriangleAroundVertexRange(std::size_t vertex_index) const;
    /**
     * Returns the starting iterator of a range over every triangle incident to
     * the vertex that has the given index within the mesh. This iterator
     * iterates on object of type Triangle.
     * @see getTriangleAroundVertexEnd
     * @see getTriangleAroundVertexRange
     */
    TriangleIterator getTriangleAroundVertexBegin(
      std::size_t vertex_index) const;
    /**
     * Returns the after end iterator of a range over every triangle incident to
     * the vertex that has the given index within the mesh.
     * @see getTriangleAroundVertexBegin
     * @see getTriangleAroundVertexRange
     */
    TriangleIterator getTriangleAroundVertexEnd(std::size_t vertex_index) const;

    /**
     * Return the length of the given edge.
     */
    TinyScalar getEdgeLength(const Edge& edge) const;

    /**
     * Returns a range over every edge of the mesh. The returned range
     * iterates over object of type Edge.
     * @getEdgeBegin
     * @getEdgeEnd
     */
    EdgeRange getEdgeRange() const;
    /**
     * Returns the starting iterator of a range over every edge in the
     * mesh. This iterator iterates on object of type Edge.
     * @see getEdgeRange
     * @see getEdgeEnd
     */
    EdgeIterator getEdgeBegin() const;
    /**
     * Returns the after end iterator of a range over every edge in the
     * mesh.
     * @see getEdgeRange
     * @see getEdgeBegin
     */
    EdgeIterator getEdgeEnd() const;

    /**
     * Returns the starting iterator of a range over the indices of the vertices
     * that are not on the a border of the mesh. A vertex is on a border if one
     * of its incident edge doesn't have two incident triangle.
     * @see getInteriorVerticesEnd
     * @see isVertexOnBorder
     * @see isInteriorVertex
     */
    InteriorVerticesIterator getInteriorVerticesBegin() const;
    /**
     * Returns the after end iterator of a range over the indices of the vertices
     * that are not on the a border of the mesh.
     * @see getInteriorVerticesBegin
     * @see isVertexOnBorder
     * @see isInteriorVertex
     */
    InteriorVerticesIterator getInteriorVerticesEnd() const;
    /**
     * Returns a range over the indices of every vertex adjacent to the vertex
     * whose index within the mesh is the given index.
     */
    VerticesAroundVertexRange getVerticesAroundVertexRange(
      std::size_t vertex_index) const;
    /**
     * Returns a range over the indices of every vertex of the mesh.
     */
    VerticesRange getVerticesRange() const;
    /**
     * Returns the position of every vertex within the mesh.
     */
    std::vector<Eigen::Matrix<TinyScalar, 3, 1>> getVertices() const noexcept;

    /**
     * Return the position of the vertex whose index within the mesh is the
     * given index.
     */
    Vertex getVertex(std::size_t index) const;
    /**
     * Returns true if and only if the vertex with the given index within the
     * mesh is on the border of the mesh. A vertex is on the border if one of
     * its incident edge has only one incident triangle. 
     * @see isInteriorVertex
     * @see getInteriorVerticesEnd
     * @see getInteriorVerticesBegin
     */
    bool isVertexOnBorder(std::size_t vertex_index) const;
    /**
     * Returns true if and only if the vertex with the given index within the
     * mesh is not on the border of the mesh. A vertex is on the border if one
     * of its incident edge has only one incident triangle.
     * @see isVertexOnBorder
     * @see getInteriorVerticesEnd
     * @see getInteriorVerticesBegin
     */
    bool isInteriorVertex(std::size_t vertex_index) const;

    /**
     * Returns the position.
     *
     * The position for a mesh is simply the concatenation of
     * the position of each vertex into one row vector. The order is
     * preserved.
     */
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& getGeneralizedPositions() const;

    //TODO: Change the name to setPositions
    /**
     * Change the position.
     *
     * It copies the given vector.
     * \see getGeneralizedPosition
     */
    void setGeneralizedPositions(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& q);

    /**
     * Change the position.
     *
     * It moves the given vector.
     * \see getGeneralizedPosition
     */
    void setGeneralizedPositions(Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>&& q);

    /**
     * Sets the position of the mesh.
     * @param positions The position of each vertices.
     */
    void setPositions(const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& positions);
    /**
     * Sets the position of the mesh.
     * @param positions The position of each vertices.
     */
    void setPositions(std::vector<Eigen::Matrix<TinyScalar, 3, 1>>&& positions);

    CollisionMesh& operator=(CollisionMesh&& other) = default;
    CollisionMesh& operator=(const CollisionMesh& other) = default;

    void writeObj(std::ofstream &of);
    std::size_t getMaterialIdentifier() const;

    std::size_t m_material_identifier;

    /**
     * Returns a function which takes the index of an edge and return the
     * associated Edge object.
     */
    std::function<Edge(const typename SurfaceMesh::Edge_index&)>
    getEdgeFromEdgeIndexFunctionObject() const {
        return std::bind(&CollisionMesh::getEdgeFromEdgeIndex, this, std::placeholders::_1);
    }

    /**
     * Returns the Edge object associated to the given edge index.
     */
    Edge getEdgeFromEdgeIndex(const typename SurfaceMesh::Edge_index& edge_index) const;

    /**
     * Returns the Triangle object associated to the given triangle index.
     */
    Triangle getTriangleFromFaceIndex(
      const typename SurfaceMesh::Face_index& face_index) const;

    /**
     * Returns a function which takes the index of a triangle and returns the
     * associated Triangle object.
     */
    std::function<Triangle(const typename SurfaceMesh::Face_index&)>
    getTriangleFromFaceIndexFunctionObject() const {
        return std::bind(&CollisionMesh::getTriangleFromFaceIndex, this, std::placeholders::_1);
    }

    // get mass
    boost::iterator_range<
      boost::transform_iterator<std::function<TinyScalar(const Triangle&)>,
                                TriangleIterator> >
      getTriangleAroundVertexMassRange(std::size_t vertex_index) const
    {
      return boost::make_iterator_range(
        boost::make_transform_iterator(getTriangleAroundVertexBegin(vertex_index),
                                       getTriangleMassFunctionObject()),
        boost::make_transform_iterator(getTriangleAroundVertexEnd(vertex_index),
                                       getTriangleMassFunctionObject()));
    }

    TinyScalar getTriangleMass(const Triangle& triangle) const
    {
        // TODO physics parameters
        // return m_physic_parameters.area_density * getTriangleArea(triangle);
        return 1.5 * getTriangleArea(triangle);
    }

    // std::function<TinyScalar(const Triangle&)> getTriangleMassFunctionObject() const;

    // std::function<TinyScalar(const Triangle&)> getTriangleMassFunctionObject() const
    // {
    //   return std::bind(&getTriangleMass, this, std::placeholders::_1);
    // }

    std::function<TinyScalar(const Triangle&)> getTriangleMassFunctionObject() const
    {
      return std::bind(&CollisionMesh::getTriangleMass, this, std::placeholders::_1);
    }


    /**
     * Returns true if the index is a null face. A null face is used in CGAL to
     * represent the absence of a triangle. For example, an edge on the border
     * of the mesh is an edge that has one null face within its incident faces.
     */
    bool isNotNullFace(const typename SurfaceMesh::Face_index& face_index) const;

    void initializeMasses();
    /**
     * A CGAL mesh that stores the topology of the mesh. The position of the
     * vertices are not stored in this data member as it would be too costly.
     * To update the position within this data member, use updateUnderlyingMesh.
     */
    SurfaceMesh m_mesh;
    /**
     * A mapping from the vertices index passed to this class methods to the
     * vertices indices within the CGAL surface mesh.  Since the CGAL surface
     * mesh does not necessarily keep the vertices in the order in which they
     * were added to the mesh, we keep a mapping.
     */
    std::vector<std::size_t> m_vertices_indices;
    /**
     * The concatenated positions of the vertices.
     */
    void SetVertexNormal();
    void initializeClothConstraints(
        std::vector<FEM::Constraint<TinyScalar, TinyConstants>*>& constraints);

    FEM::Constraint<TinyScalar, TinyConstants>* getEdgeConstraint(const Edge& edge) const;

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> m_positions;
    std::vector<int> m_triangles;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> m_speed;
    std::vector<TinyScalar> m_masses;
    ClothPhysicParameters<TinyScalar, TinyConstants> mParameter;
    std::vector<Eigen::Vector3i>            mContours;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                         mVertexNormal;
};

// Include implementation of templated methods.

template <typename TinyScalar, typename TinyConstants> 
CollisionMesh<TinyScalar, TinyConstants>::
CollisionMesh(const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& positions,
           const std::vector<int>& triangles)
  : m_vertices_indices(positions.size()), m_positions(positions.size() * 3),
  m_triangles(triangles)
{
    using VertexIndex = typename CGAL::Surface_mesh<Eigen::Matrix<TinyScalar, 3, 1> >::Vertex_index;

    VertexIndex current_vertex_index;
    for (std::size_t i = 0; i < positions.size(); ++i)
    {
        current_vertex_index = m_mesh.add_vertex(positions[i]);
        // We want to be sure that the indices will stay inside the generalized
        // position ranged.
        assert(current_vertex_index < positions.size());
        m_vertices_indices[i] = current_vertex_index;
        getVector3dBlock<TinyScalar, TinyConstants>(this->m_positions, current_vertex_index) = positions[i];
    }
    std::cout << triangles.size() << " triangles.size()\n";
    std::cout << positions.size() << " positions.size()\n";

    for (std::size_t i = 0; i < triangles.size() / 3; ++i)
    {
        // TODO: Auxiliary function to get triangle edges
        m_mesh.add_face(
          (typename SurfaceMesh::Vertex_index)m_vertices_indices[triangles[3 * i]],
          (typename SurfaceMesh::Vertex_index)m_vertices_indices[triangles[3 * i + 1]],
          (typename SurfaceMesh::Vertex_index)m_vertices_indices[triangles[3 * i + 2]]);
    }

    for(int i=0;i<triangles.size() / 3;i++)
    {
        Eigen::Vector3i t_1 = Eigen::Vector3i(
            triangles[3 * i], triangles[3 * i + 1], triangles[3 * i + 2]);
        mContours.push_back(t_1);
        Eigen::Vector3i t_2 = Eigen::Vector3i(
            triangles[3 * i], triangles[3 * i + 2], triangles[3 * i + 1]);
        mContours.push_back(t_2);
    }

    initializeMasses();
}

template <typename TinyScalar, typename TinyConstants> 
void CollisionMesh<TinyScalar, TinyConstants>::
initializeMasses()
{
  m_masses.resize(getNumberOfVertices());
  for (std::size_t vId = 0; vId < m_masses.size(); ++vId)
  {
    auto triangle_around_vertex_mass_range =
      getTriangleAroundVertexMassRange(vId);
    m_masses[vId] = std::accumulate(boost::begin(triangle_around_vertex_mass_range),
                            boost::end(triangle_around_vertex_mass_range),
                            TinyConstants::zero()) / 3;
  }
}


template <typename TinyScalar, typename TinyConstants> 
CollisionMesh<TinyScalar, TinyConstants>::
CollisionMesh(const SurfaceMesh& mesh) : m_mesh(mesh)
{
    m_vertices_indices.resize(m_mesh.number_of_vertices());
    std::iota(m_vertices_indices.begin(), m_vertices_indices.end(), 0);
    for (const auto& vertex_index : m_mesh.vertices())
    {
        getVector3dBlock<TinyScalar, TinyConstants>(m_positions, vertex_index) = m_mesh.point(vertex_index);
    }
}

template <typename TinyScalar, typename TinyConstants> 
const typename CollisionMesh<TinyScalar, TinyConstants>::SurfaceMesh&
CollisionMesh<TinyScalar, TinyConstants>::
getUnderlyingMesh() const
{
    return m_mesh;
}

template <typename TinyScalar, typename TinyConstants> 
std::size_t CollisionMesh<TinyScalar, TinyConstants>::
getMaterialIdentifier() const
{
    return m_material_identifier;
}


template <typename TinyScalar, typename TinyConstants> 
void CollisionMesh<TinyScalar, TinyConstants>::
updateUnderlyingMesh()
{
    for (auto& vertex_index : m_mesh.vertices())
    {
        m_mesh.point(vertex_index) = getVector3dBlock<TinyScalar, TinyConstants>(m_positions, vertex_index);
    }
}

template <typename TinyScalar, typename TinyConstants> 
std::size_t CollisionMesh<TinyScalar, TinyConstants>::
getNumberOfVertices() const
{
    return getVector3dSize<TinyScalar, TinyConstants>(getGeneralizedPositions());
}

template <typename TinyScalar, typename TinyConstants> 
std::size_t CollisionMesh<TinyScalar, TinyConstants>::
getNumberOfTriangles() const noexcept
{
    return m_mesh.number_of_faces();
}

template <typename TinyScalar, typename TinyConstants> 
void CollisionMesh<TinyScalar, TinyConstants>::
writeObj(std::ofstream &of)
{
  for (size_t vId = 0u; vId < getNumberOfVertices(); ++vId)
  {
    const Eigen::Matrix<TinyScalar, 3, 1> v = getVertex(vId);
    of << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
  }
  for (size_t fId = 0u; fId < getNumberOfTriangles(); ++fId)
  {
    const CollisionMesh::Triangle t = getTriangle(fId);
    of << "f " << t.vertex_indices[0] + 1
       << " " << t.vertex_indices[1] + 1
       << " " << t.vertex_indices[2] + 1 << std::endl;
  }
}


template <typename TinyScalar, typename TinyConstants> 
TinyScalar CollisionMesh<TinyScalar, TinyConstants>::
getTriangleArea(const Triangle& triangle) const
{
    Eigen::Matrix<double, 3, 1> v1 = helper::to_eigen_double<TinyScalar, TinyConstants>(
        getVertex(triangle.vertex_indices[0]));
    Eigen::Matrix<double, 3, 1> v2 = helper::to_eigen_double<TinyScalar, TinyConstants>(
        getVertex(triangle.vertex_indices[1]));
    Eigen::Matrix<double, 3, 1> v3 = helper::to_eigen_double<TinyScalar, TinyConstants>(
        getVertex(triangle.vertex_indices[2]));
    double tmp = CGAL::squared_area<EigenKernel<double, DoubleUtils> >(v1, v2, v3);
    // TinyScalar tmp = CGAL::squared_area<EigenKernel<TinyScalar, TinyConstants> >(
    //     getVertex(triangle.vertex_indices[0]),
    //     getVertex(triangle.vertex_indices[1]),
    //     getVertex(triangle.vertex_indices[2]));
    // TinyScalar tmp = TinyConstants::zero();
    // TODO conversion CGAL::Null_tag
    return TinyConstants::sqrt1(TinyConstants::scalar_from_double(tmp));
}

template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, 3, 1> CollisionMesh<TinyScalar, TinyConstants>::
getTriangleNormal(const Triangle& triangle) const
{
    return (getVertex(triangle.vertex_indices[1]) - getVertex(triangle.vertex_indices[0]))
      .cross(getVertex(triangle.vertex_indices[2]) - getVertex(triangle.vertex_indices[0]))
      .normalized();
}

template <typename TinyScalar, typename TinyConstants> 
CGAL::Bbox_3 CollisionMesh<TinyScalar, TinyConstants>::
getTriangleBoundingBox(const Triangle& triangle) const
{
    // return CGAL::Bbox_3();
    // TODO conversion CGAL::Null_tag


    Eigen::Matrix<double, 3, 1> v1 = helper::to_eigen_double<TinyScalar, TinyConstants>(
        getVertex(triangle.vertex_indices[0]));
    Eigen::Matrix<double, 3, 1> v2 = helper::to_eigen_double<TinyScalar, TinyConstants>(
        getVertex(triangle.vertex_indices[1]));
    Eigen::Matrix<double, 3, 1> v3 = helper::to_eigen_double<TinyScalar, TinyConstants>(
        getVertex(triangle.vertex_indices[2]));

    return typename EigenKernel<double, DoubleUtils>::Triangle_3(v1, v2, v3).bbox();
}

/**
 *  @brief Return the position of the triangle vertices in the column of a matrix.
 */
template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, 3, 3>
CollisionMesh<TinyScalar, TinyConstants>::
getTriangle(const Triangle& triangle) const noexcept
{
    Eigen::Matrix<TinyScalar, 3, 3> result;
    result.col(0) = getVertex(triangle.vertex_indices[0]);
    result.col(1) = getVertex(triangle.vertex_indices[1]);
    result.col(2) = getVertex(triangle.vertex_indices[2]);
    return result;
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::Triangle
CollisionMesh<TinyScalar, TinyConstants>::
getTriangle(std::size_t triangle_index) const noexcept
{
    return getTriangleFromFaceIndex(typename SurfaceMesh::Face_index(triangle_index));
};

template <typename TinyScalar, typename TinyConstants> 
TinyScalar CollisionMesh<TinyScalar, TinyConstants>::
getEdgeLength(const Edge& edge) const
{
    return (getVertex(edge.vertex1_index) - getVertex(edge.vertex2_index)).norm();
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::AllTriangleRange
CollisionMesh<TinyScalar, TinyConstants>::
getTriangleRange() const noexcept
{
    return AllTriangleRange(getTriangleBegin(), getTriangleEnd());
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::AllTriangleIterator
CollisionMesh<TinyScalar, TinyConstants>::
getTriangleBegin() const noexcept
{
    return AllTriangleIterator(m_mesh.faces_begin(), getTriangleFromFaceIndexFunctionObject());
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::AllTriangleIterator
CollisionMesh<TinyScalar, TinyConstants>::
getTriangleEnd() const noexcept
{
    return AllTriangleIterator(m_mesh.faces_end(), getTriangleFromFaceIndexFunctionObject());
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::TriangleRange
CollisionMesh<TinyScalar, TinyConstants>::
getTriangleAroundVertexRange(std::size_t vertex_index) const
{
    return TriangleRange(getTriangleAroundVertexBegin(vertex_index),
                         getTriangleAroundVertexEnd(vertex_index));
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::TriangleIterator
CollisionMesh<TinyScalar, TinyConstants>::
getTriangleAroundVertexBegin(std::size_t vertex_index) const
{
    using FaceAroundTargetExceptSomeIterator =
      boost::filter_iterator<std::function<bool(const typename SurfaceMesh::Face_index&)>,
                             typename SurfaceMesh::Face_around_target_iterator>;

    auto triangle_around_vertex_range =
      m_mesh.faces_around_target(m_mesh.halfedge((
        typename SurfaceMesh::Vertex_index)vertex_index));
    return TriangleIterator(FaceAroundTargetExceptSomeIterator(
                              std::bind(&CollisionMesh::isNotNullFace, this, std::placeholders::_1),
                              boost::begin(triangle_around_vertex_range),
                              boost::end(triangle_around_vertex_range)),
                            getTriangleFromFaceIndexFunctionObject());
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::TriangleIterator
CollisionMesh<TinyScalar, TinyConstants>::
getTriangleAroundVertexEnd(std::size_t vertex_index) const
{
    using FaceAroundTargetExceptSomeIterator =
      boost::filter_iterator<std::function<bool(const typename SurfaceMesh::Face_index&)>,
                             typename SurfaceMesh::Face_around_target_iterator>;

    auto triangle_around_vertex_range =
      m_mesh.faces_around_target(m_mesh.halfedge((typename SurfaceMesh::Vertex_index)vertex_index));
    return TriangleIterator(FaceAroundTargetExceptSomeIterator(
                              std::bind(&CollisionMesh::isNotNullFace, this, std::placeholders::_1),
                              boost::end(triangle_around_vertex_range),
                              boost::end(triangle_around_vertex_range)),
                            getTriangleFromFaceIndexFunctionObject());
}

template <typename TinyScalar, typename TinyConstants>
typename CollisionMesh<TinyScalar, TinyConstants>::EdgeRange 
CollisionMesh<TinyScalar, TinyConstants>::
getEdgeRange() const
{
    return EdgeRange(getEdgeBegin(), getEdgeEnd());
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::EdgeIterator
CollisionMesh<TinyScalar, TinyConstants>::
getEdgeBegin() const
{
    return EdgeIterator(m_mesh.edges_begin(), getEdgeFromEdgeIndexFunctionObject());
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::EdgeIterator
CollisionMesh<TinyScalar, TinyConstants>::
getEdgeEnd() const
{
    return EdgeIterator(m_mesh.edges_end(), getEdgeFromEdgeIndexFunctionObject());
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::InteriorVerticesIterator
CollisionMesh<TinyScalar, TinyConstants>::
getInteriorVerticesBegin() const
{
    return InteriorVerticesIterator(std::bind(&CollisionMesh::isInteriorVertex, this, 
                                    std::placeholders::_1),
                                    boost::begin(m_mesh.vertices()),
                                    boost::end(m_mesh.vertices()));
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::InteriorVerticesIterator 
CollisionMesh<TinyScalar, TinyConstants>::
getInteriorVerticesEnd() const
{
    return InteriorVerticesIterator(std::bind(&CollisionMesh::isInteriorVertex, this, std::placeholders::_1),
                                    boost::end(m_mesh.vertices()),
                                    boost::end(m_mesh.vertices()));
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::VerticesAroundVertexRange 
CollisionMesh<TinyScalar, TinyConstants>::
getVerticesAroundVertexRange(std::size_t vertex_index) const
{
    return m_mesh.vertices_around_target(m_mesh.halfedge((typename SurfaceMesh::Vertex_index)vertex_index));
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::VerticesRange 
CollisionMesh<TinyScalar, TinyConstants>::
getVerticesRange() const
{
    return m_mesh.vertices();
}

template <typename TinyScalar, typename TinyConstants> 
std::vector<Eigen::Matrix<TinyScalar, 3, 1>> CollisionMesh<TinyScalar, TinyConstants>::
getVertices() const noexcept
{
    std::vector<Eigen::Matrix<TinyScalar, 3, 1>> result(getNumberOfVertices());
    for (std::size_t vertex_index : getVerticesRange())
    {
        result[vertex_index] = getVertex(vertex_index);
    }
    return result;
}

template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::Vertex 
CollisionMesh<TinyScalar, TinyConstants>::
getVertex(std::size_t index) const
{
    // printf("getVertex %d %d\n", index, getGeneralizedPositions().size());
    return getVector3dBlock<TinyScalar, TinyConstants>(getGeneralizedPositions(), index);
}

template <typename TinyScalar, typename TinyConstants> 
bool CollisionMesh<TinyScalar, TinyConstants>::
isVertexOnBorder(
    std::size_t vertex_index) const
{
    return m_mesh.is_border((typename SurfaceMesh::Vertex_index)vertex_index);
}

template <typename TinyScalar, typename TinyConstants> 
bool CollisionMesh<TinyScalar, TinyConstants>::
isInteriorVertex(std::size_t vertex_index) const
{
    return !isVertexOnBorder(vertex_index);
}

template <typename TinyScalar, typename TinyConstants> 
const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& CollisionMesh<TinyScalar, TinyConstants>::
getGeneralizedPositions() const
{
    return m_positions;
}

template <typename TinyScalar, typename TinyConstants> 
void CollisionMesh<TinyScalar, TinyConstants>::
setGeneralizedPositions(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& q)
{
    m_positions = q;
}

template <typename TinyScalar, typename TinyConstants> 
void CollisionMesh<TinyScalar, TinyConstants>::
setGeneralizedPositions(Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>&& q)
{
    m_positions = std::move(q);
}

template <typename TinyScalar, typename TinyConstants> 
void CollisionMesh<TinyScalar, TinyConstants>::
setPositions(const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& positions)
{
    setPositions(std::vector<Eigen::Matrix<TinyScalar, 3, 1>>(positions));
}

template <typename TinyScalar, typename TinyConstants> 
void CollisionMesh<TinyScalar, TinyConstants>::
setPositions(std::vector<Eigen::Matrix<TinyScalar, 3, 1>>&& positions)
{
    assert(positions.size() == getNumberOfVertices());
    for (std::size_t i = 0; i < positions.size(); ++i)
    {
        getVector3dBlock<TinyScalar, TinyConstants>(m_positions, m_vertices_indices[i]) = positions[i];
    }
}


template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::Edge 
CollisionMesh<TinyScalar, TinyConstants>::
getEdgeFromEdgeIndex(const typename SurfaceMesh::Edge_index& edge_index) const
{
    return Edge{ m_mesh.vertex(edge_index, 0), m_mesh.vertex(edge_index, 1) };
}


template <typename TinyScalar, typename TinyConstants> 
typename CollisionMesh<TinyScalar, TinyConstants>::Triangle 
CollisionMesh<TinyScalar, TinyConstants>::
getTriangleFromFaceIndex(const typename SurfaceMesh::Face_index& face_index) const
{
    auto vertices_around_face_range = m_mesh.vertices_around_face(m_mesh.halfedge(face_index));

    Triangle triangle;
    std::copy(boost::begin(vertices_around_face_range),
              boost::end(vertices_around_face_range),
              triangle.vertex_indices.begin());

    return triangle;
}

template <typename TinyScalar, typename TinyConstants> 
bool CollisionMesh<TinyScalar, TinyConstants>::
isNotNullFace(const typename SurfaceMesh::Face_index& face_index) const
{
    return face_index != m_mesh.null_face();
}

template <typename TinyScalar, typename TinyConstants> 
void CollisionMesh<TinyScalar, TinyConstants>::
SetVertexNormal()
{

    mVertexNormal.resize(m_positions.size());
    mVertexNormal.setZero();
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> cnt_list(m_positions.size()/3);
    cnt_list.setZero();

    for(int i=0;i<mContours.size();i++)
    {
        int i0 = mContours[i][0];
        int i1 = mContours[i][1];
        int i2 = mContours[i][2];

        Eigen::Matrix<TinyScalar, 3, 1> p0 = getVector3dBlock<TinyScalar, TinyConstants>(getGeneralizedPositions(), i0);
        Eigen::Matrix<TinyScalar, 3, 1> p1 = getVector3dBlock<TinyScalar, TinyConstants>(getGeneralizedPositions(), i1);
        Eigen::Matrix<TinyScalar, 3, 1> p2 = getVector3dBlock<TinyScalar, TinyConstants>(getGeneralizedPositions(), i2);

        Eigen::Matrix<TinyScalar, 3, 1> face_normal;
        face_normal = (p1-p0).cross(p2-p1);

        mVertexNormal.template block<3,1>(3*i0,0) = cnt_list[i0]/(cnt_list[i0]+1)*mVertexNormal.template block<3,1>(3*i0,0) + 1/(cnt_list[i0]+1)*face_normal;
        mVertexNormal.template block<3,1>(3*i1,0) = cnt_list[i1]/(cnt_list[i1]+1)*mVertexNormal.template block<3,1>(3*i1,0) + 1/(cnt_list[i1]+1)*face_normal;
        mVertexNormal.template block<3,1>(3*i2,0) = cnt_list[i2]/(cnt_list[i2]+1)*mVertexNormal.template block<3,1>(3*i2,0) + 1/(cnt_list[i2]+1)*face_normal;

        cnt_list[i0] += 1;
        cnt_list[i1] += 1;
        cnt_list[i2] += 1;
    }
}


template <typename TinyScalar, typename TinyConstants> 
void CollisionMesh<TinyScalar, TinyConstants>::
initializeClothConstraints(
    std::vector<FEM::Constraint<TinyScalar, TinyConstants>*>& constraints)
{
  std::cout << " initializeClothConstraints\n";
  std::transform(
    getEdgeBegin(),
    getEdgeEnd(),
    std::back_inserter(constraints),
    std::bind(&CollisionMesh::getEdgeConstraint, this, std::placeholders::_1));
  // std::transform(
  //   getInteriorVerticesBegin(),
  //   getInteriorVerticesEnd(),
  //   std::back_inserter(constraints),
  //   std::bind(&CollisionMesh::getVertexConstraint, this, std::placeholders::_1));
}

// template <typename TinyScalar, typename TinyConstants> 
// std::unique_ptr<FEM::Constraint<TinyScalar, TinyConstants> > 
// CollisionMesh<TinyScalar, TinyConstants>::
// getVertexConstraint(std::size_t vertex_index) const
// {
//   auto vertices_around_vertex = getVerticesAroundVertexRange(vertex_index);
//   std::vector<std::size_t> vertex_indices;
//   vertex_indices.push_back(vertex_index);
//   vertex_indices.insert(vertex_indices.end(),
//                         boost::begin(vertices_around_vertex),
//                         boost::end(vertices_around_vertex));

//   return std::make_unique<BendingConstraint>(
//     vertex_indices,
//     getVertex(vertex_index),
//     std::vector<Eigen::Matrix<TinyScalar, 3, 1>>(
//       boost::make_transform_iterator(
//         boost::begin(vertices_around_vertex),
//         std::bind(&Mesh::getVertex, this, std::placeholders::_1)),
//       boost::make_transform_iterator(
//         boost::end(vertices_around_vertex),
//         std::bind(&Mesh::getVertex, this, std::placeholders::_1))),
//     m_physic_parameters.bend);
//     // m_physic_parameters.bend);
// }


template <typename TinyScalar, typename TinyConstants> 
FEM::Constraint<TinyScalar, TinyConstants>*
CollisionMesh<TinyScalar, TinyConstants>::
getEdgeConstraint(const Edge& edge) const
{
  TinyScalar rest_length = getEdgeLength(edge);
  SpringConstraint<TinyScalar, TinyConstants>* c = 
    new SpringConstraint<TinyScalar, TinyConstants>(
    edge.vertex1_index, edge.vertex2_index, rest_length, 200);
  return c;
    // edge.vertex1_index, edge.vertex2_index, rest_length, m_physic_parameters.stretch);
}


};
#endif
