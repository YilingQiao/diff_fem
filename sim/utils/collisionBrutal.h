#ifndef __COLLISION_BRUTAL__H__
#define __COLLISION_BRUTAL__H__
#include <CGAL/Surface_mesh.h>
#include <boost/iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/sub_range.hpp>
#include <Eigen/Core>
#include <functional>
#include <tiny_obj_loader.h>
#include <vector>
#include <fstream>
#include <CGAL/box_intersection_d.h>

#include "arcsim/collision.hpp"
#include "../fem/Mesh/CollisionMesh.h"
#include "timer.h"
namespace FEM
{ 

/// @brief Struct containing the data needed to handle a collision
template <typename TinyScalar, typename TinyConstants> 
struct CollisionInfo 
{
  /// @brief (Generalized) index of the vertex in contact.
  size_t vertex_index;
  /// @brief Position of the contact point.
  Eigen::Matrix<TinyScalar, 3, 1> contact_point;
  /// @brief Normal at the contact point.
  Eigen::Matrix<TinyScalar, 3, 1> normal;
  /// Speed of the contact point.
  Eigen::Matrix<TinyScalar, 3, 1> speed;
  /// @brief Friction coefficient of the two material in collision.
  TinyScalar friction_coefficient;
  CollisionInfo(
      int vi, 
      const Eigen::Matrix<TinyScalar, 3, 1> &cp,
      const Eigen::Matrix<TinyScalar, 3, 1> &n, 
      const Eigen::Matrix<TinyScalar, 3, 1> &s, TinyScalar fc) {
    vertex_index = vi;
    contact_point = cp;
    normal = n;
    speed = s;
    friction_coefficient = fc;
  }
  CollisionInfo(){}
  static const TinyScalar default_friction;
};

template<typename TinyScalar, typename TinyConstants>
const TinyScalar CollisionInfo<TinyScalar, TinyConstants>::default_friction = 1;

template <typename TinyScalar, typename TinyConstants> 
struct SelfForceToAdd
{
    size_t id_plus;
    std::array<int,3> id_minus;
    Eigen::Matrix<TinyScalar, 3, 1> alpha;
    // size_t id_minus;
    Eigen::Matrix<TinyScalar, 3, 1> force;
};

/**
 * Stores the date required to handle a self-collision. Currently,
 * self-collision is only node-node.
 */
template <typename TinyScalar, typename TinyConstants> 
struct SelfCollisionInfo : public CollisionInfo<TinyScalar, TinyConstants>
{
  /**
   * Index of the vertices that make up the triangle with which the vertex
   * collided.
   */
  std::array<int, 3> face_indices;
  /**
   * The barycentric coordinate of the point in the face with which the vertex
   * has collided. Currently, since we only handle node-node self-collision,
   * this vector has a coordiante set to 1 and the other to 0.
   */
  Eigen::Matrix<TinyScalar, 3, 1> barycentric_coordinates;
  size_t vertex_index_b;

  std::vector<size_t> vertex_index_list;
};


template <typename TinyScalar, typename TinyConstants> 
struct IndexBox
{
    CGAL::Box_intersection_d::Box_d<double, 3u> m_box;
    std::size_t m_index;
    using NT = double;
    using ID = std::size_t; 
    
    IndexBox(const IndexBox& other) = default;
    IndexBox(const CGAL::Bbox_3& box, std::size_t index) : m_box(box), m_index(index) {} 
    static int dimension() { return 3u; };
    ID id() const noexcept { return m_box.id(); }
    NT min_coord(int d) const noexcept { return m_box.min_coord(d); }
    NT max_coord(int d) const noexcept { return m_box.max_coord(d); }
    IndexBox& operator=(const IndexBox& other) = default;
};


template <typename TinyScalar, typename TinyConstants> 
class CollisionBrutal
{
    using FrictionCoefficientTable = std::vector<std::vector<TinyScalar>>;
public:
  CollisionBrutal<TinyScalar, TinyConstants>(
    CollisionMesh<TinyScalar, TinyConstants>* collisionMesh, bool handle_self_collision);

  
  bool BaseUpdateCollisions(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& X);
  void ClearCollisions();
  void LocalCollision() noexcept;
  std::vector<CollisionInfo<TinyScalar, TinyConstants> > getCollisionsInfo(
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& positions) const;
  std::vector<SelfCollisionInfo<TinyScalar, TinyConstants> > getSelfCollisionsInfo(
    TinyScalar tol2) ;
  std::vector<SelfCollisionInfo<TinyScalar, TinyConstants> > checkMeshSelfCollisions(
    const CollisionMesh<TinyScalar, TinyConstants>* mesh, TinyScalar tolerance) noexcept;
  void computeBasis(const Eigen::Matrix<TinyScalar, 3, 1>& normal, Eigen::Matrix<TinyScalar, 3, 3>& basis);
  void computeCollisionComputationOrder() noexcept;
  void computeSelfCollisionGraph() noexcept;
  void fillCollisionComputationOrder(
    const std::vector<std::size_t>& starting_vertices,
    std::vector<bool>& visited,
    std::size_t initial_computation_index) noexcept;
  std::size_t getSelfCollisionContactPointIndex(
    const SelfCollisionInfo<TinyScalar, TinyConstants>& self_collision_info) const noexcept;
  bool vertexIsVisited(std::size_t vertex_index, const std::vector<bool>& visited) const noexcept;
  std::vector<std::vector<
  // std::pair<typename CollisionMesh<TinyScalar, TinyConstants>::Triangle, Eigen::Matrix<TinyScalar, 3, 1>> > > 
  std::pair<Impact*, Eigen::Matrix<TinyScalar, 3, 1>> > > 
  checkMeshAllCollision(
    const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& vertices, 
    const CollisionMesh<TinyScalar, TinyConstants>* mesh, 
    TinyScalar tolerance);
  std::vector<std::vector<
    typename CollisionMesh<TinyScalar, TinyConstants>::Triangle> > 
    getMeshPotentialTriangleCollision(
    const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& vertices,
    const CollisionMesh<TinyScalar, TinyConstants>* mesh,
    TinyScalar tolerance);
  CGAL::Bbox_3 getToleranceBoundingBox(const Eigen::Matrix<TinyScalar, 3, 1>& vertex, TinyScalar tolerance);
  void closestPointToTriangle(const Eigen::Matrix<TinyScalar, 3, 1>& p,
    const Eigen::Matrix<TinyScalar, 3, 3>& triangle,
    TinyScalar tol2,
    Eigen::Matrix<TinyScalar, 3, 1>& closest_point,
    Eigen::Matrix<TinyScalar, 3, 1>* barycentric_coordinates_ptr) noexcept;
  std::size_t isTriangleVertexIndexIn(
    const typename CollisionMesh<TinyScalar, TinyConstants>::Triangle& triangle, 
    const std::vector<size_t>& vertices_indices);


  std::vector<int> m_collision_numbers;
  // TinyScalar m_self_collision_tol2 = 1e-6; 
  // TinyScalar m_self_collision_tol2 = 0.0025;
  double m_self_collision_tol2 = 1e-4;    
  std::vector<TinyScalar> m_step_times;
  std::vector<TinyScalar> m_rhs_times;
  std::vector<TinyScalar> m_global_times;
  std::vector<TinyScalar> m_iteration_times;
  std::vector<TinyScalar> m_collision_detection_times;
  std::vector<TinyScalar> m_self_collision_detection_times;
  std::vector<size_t> m_self_collision_numbers;
  CollisionMesh<TinyScalar, TinyConstants>* mCollisionMesh;
  
  
  FrictionCoefficientTable m_friction_coefficients;
  //@brief Collision info should be kept sorted in respect to the vertex indices
  //TODO: use a better data structure for sorted collections.
  std::vector<CollisionInfo<TinyScalar, TinyConstants>> m_collisions_infos;
  std::vector<SelfCollisionInfo<TinyScalar, TinyConstants>> m_self_collisions_infos;
  TinyScalar m_damping_coefficient;
  /// @brief
  Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> m_damping;
    bool m_handle_self_collision = true;

  // frecitionestimation

    std::vector<Eigen::Matrix<TinyScalar, 3, 3>> m_local_contact_basis;

    /**
     * Represent an edge in the self-collision graph adjacency list.
     * @see m_self_collision_graph
     */
    struct SelfCollisionGraphNeighboor
    {
        /**
         * The index of vertex adjacent through this edge.
         */
        std::size_t index;
        /**
         * The index of the self-collision represented by this index.
         */
        std::size_t collision_index;
    };
    /**
     * An adjacency list representing the self-collision graph. The vertices of
     * this graph are the vertices of the mesh. The edges are the self-collision
     * between the vertices.
     * @see SelfCollisionGraphNeighboor
     * @see computeSelfCollisionGraph
     */
    std::vector<std::vector<SelfCollisionGraphNeighboor>> m_self_collision_graph;
    std::vector<Eigen::Matrix<TinyScalar, 3, 3>> m_local_self_contact_basis;
    /**
     * The order in which the computation of repulsive force should be made.
     * The first element of this vector contains the indices of the
     * self-collision whose repulsion force should be computed first, the second
     * element those whose repulsion force that should be computed in second. So
     * on and so forth. Its value can be computed through
     * computeSelfCollisionGraph. The self-collision on a same level can be
     * safelly computed in parallel.
     */
    std::vector<std::vector<std::size_t>> m_collision_computation_order;

    /// @brief (IO)
    std::vector<TinyScalar> m_rhs_base_times;
    /// @brief (IO)
    std::vector<TinyScalar> m_friction_times;
    /// @brief (IO)
    std::vector<TinyScalar> m_self_friction_times;
    std::vector<TinyScalar> m_self_collision_ordering_times;

    /// Better Row or column ?
    // column is better to write in
    // row is better to add to rhs at the end
    using ForceType = Eigen::Matrix<TinyScalar, 3, -1, Eigen::RowMajor>;

    unsigned int m_current_index;
    unsigned int m_next_index;
    /// @brief Storage needed to handle self friction
    ForceType m_contact_forces[2];
    /// @brief Storage needed to handle self friction
    ForceType m_self_contact_forces[2];
    /// @brief Storage needed to handle self friction
    ForceType m_self_contact_repercusion_forces[2];
    /// @brief Storage needed to handle self friction
    Eigen::Matrix<TinyScalar, 3, Eigen::Dynamic> m_alpha[2];

    Eigen::Matrix<TinyScalar, 3, Eigen::Dynamic> m_remember_self_contact_forces;
    /// @brief ...
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> m_generalized_positions;
    /// @brief ...
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> m_generalized_speeds;
    ArcsimMesh *mMesh, *mObsMesh;
    vector<Impact> impacts, obsImpacts;
};

template <typename TinyScalar, typename TinyConstants> 
CollisionBrutal<TinyScalar, TinyConstants>::
CollisionBrutal(CollisionMesh<TinyScalar, TinyConstants>* collisionMesh,
  bool handle_self_collision)
:mCollisionMesh(collisionMesh), m_handle_self_collision(handle_self_collision)
{

  double air_damping = 0.001;
  int nVertices = mCollisionMesh->getNumberOfVertices();
  m_damping_coefficient = TinyConstants::scalar_from_double(std::max(air_damping, 0.)) ;

  m_damping.resize(nVertices) ;
  m_self_collision_graph.resize(nVertices) ;

  // m_t_n = Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>::Zero(m_nDofs);
  // m_current_next_position = Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>::Zero(m_nDofs);
}


template <typename TinyScalar, typename TinyConstants> 
void CollisionBrutal<TinyScalar, TinyConstants>::
ClearCollisions() {
  impacts.clear();
  obsImpacts.clear();
  m_collisions_infos.clear();
  m_self_collisions_infos.clear();
}


template <typename TinyScalar, typename TinyConstants> 
bool CollisionBrutal<TinyScalar, TinyConstants>::
BaseUpdateCollisions(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& X) 
{
  TIMER_START(collision_detection);
  
  vector<Impact> newimpacts, newobsImpacts;
  int tmp = m_self_collisions_infos.size();
  boost::tie(newimpacts, newobsImpacts) = collision_detection(mMesh, mObsMesh);
  //add and unique
  for (int i = 0; i < newimpacts.size(); ++i) {
    int j = 0;
    for (j = 0; j < impacts.size(); ++j) {
      bool flag = true;
      for (int k = 0; k < 4; ++k)
        if (impacts[j].nodes[k]->index!=newimpacts[i].nodes[k]->index) {
          flag = false;
          break;
        }
      if (flag) {
        impacts[j] = newimpacts[i];
        break;
      }
    }
    if (j == impacts.size())
      impacts.push_back(newimpacts[i]);
  }
  for (int i = 0; i < newobsImpacts.size(); ++i) {
    int j = 0;
    for (j = 0; j < obsImpacts.size(); ++j) {
      bool flag = true;
      for (int k = 0; k < 4; ++k)
        if (obsImpacts[j].nodes[k]->index!=newobsImpacts[i].nodes[k]->index) {
          flag = false;
          break;
        }
      if (flag) {
        obsImpacts[j] = newobsImpacts[i];
        break;
      }
    }
    if (j == obsImpacts.size())
      obsImpacts.push_back(newobsImpacts[i]);
  }
  m_collisions_infos = getCollisionsInfo(X);
  
  if (m_collisions_infos.size() > 0) 
  {
    m_collision_numbers.push_back(m_collisions_infos.size());
  }

  m_damping.setConstant(m_damping_coefficient);
  if (m_damping_coefficient > 0.)
  {
// #pragma omp parallel
    for (size_t cId = 0u; cId < m_collisions_infos.size(); ++cId)
    {
      m_damping[m_collisions_infos[cId].vertex_index] = 0.;
    } // cId
  }
  
  const TinyScalar duration_collision_detection =
    TIMER_DURATION(collision_detection, microseconds);
  m_collision_detection_times.push_back(duration_collision_detection);
#ifdef TIMER_PRINT
  std::cout << "# Collision detection          : "
            << duration_collision_detection
            << " µs" << std::endl;
  std::cout << "Number of collision            : "
            << m_collisions_infos.size() << std::endl;
#endif // TIMER_PRINT

  if (m_handle_self_collision)
  {
    TIMER_START(self_collision_detection);
    
    m_self_collisions_infos =
      getSelfCollisionsInfo(m_self_collision_tol2);

    if (m_self_collisions_infos.size() > 0) 
    {
      m_self_collision_numbers.push_back(m_self_collisions_infos.size());
    }

    if (m_damping_coefficient > 0.)
    {
#pragma omp parallel
      for (size_t scId = 0u; scId < m_self_collisions_infos.size(); ++scId)
      {
        // std::cout << m_damping.size() << " m_damping.size() \n";
        // std::cout << m_self_collisions_infos[scId].vertex_index << 
        //   " m_self_collisions_infos[scId].vertex_index \n";
        // m_damping[m_self_collisions_infos[scId].vertex_index] = 0.;
        // // other vertex
        // const Eigen::Matrix<TinyScalar, 3, 1> &alpha = m_self_collisions_infos[scId].barycentric_coordinates;
        // const std::array<size_t, 3> &nId = m_self_collisions_infos[scId].face_indices;
        // const size_t ovId =
        //   (alpha[0] > alpha[1]) ?
        //   ((alpha[0] > alpha[2]) ? nId[0] : nId[2]) :
        //   ((alpha[1] > alpha[2]) ? nId[1] : nId[2]);
        for (int k = 0; k < m_self_collisions_infos[scId].vertex_index_list.size(); ++k)
          m_damping[m_self_collisions_infos[scId].vertex_index_list[k]] = 0.;

        // const size_t ovId = m_self_collisions_infos[scId].vertex_index_b;
        // m_damping[ovId] = 0.;
      } // vId
    }

    
    const TinyScalar duration_self_collision_detection =
      TIMER_DURATION(self_collision_detection, microseconds);
    m_self_collision_detection_times.push_back(duration_self_collision_detection);
#ifdef TIMER_PRINT
    std::cout << "# Self-collision detection     : "
              << duration_self_collision_detection
              << " µs" << std::endl;
    std::cout << "Number of self collision       : " << m_self_collisions_infos.size() << std::endl;
#endif // TIMER_PRINT

    
  }
  
  LocalCollision();
  // std::cout << tmp << " " << m_self_collisions_infos.size() << std::endl;
  return newobsImpacts.empty() && tmp == m_self_collisions_infos.size(); //newobsImpacts.empty();
}

template <typename TinyScalar, typename TinyConstants> 
void CollisionBrutal<TinyScalar, TinyConstants>::
computeBasis(const Eigen::Matrix<TinyScalar, 3, 1>& normal, Eigen::Matrix<TinyScalar, 3, 3>& basis)
{
  
  const Eigen::Matrix<TinyScalar, 3, 1> test1(1., 0., 0.);
  const Eigen::Matrix<TinyScalar, 3, 1> test2(0., 1., 0.);
    Eigen::Matrix<TinyScalar, 3, 1> tangent = normal.cross(test1);
    if (tangent.norm() < 1.e-15)
    {
        tangent = normal.cross(test2);
    }
    tangent.normalize();
    const Eigen::Matrix<TinyScalar, 3, 1> bitangent = normal.cross(tangent).normalized();

    basis.col(0) = normal.normalized();
    basis.col(1) = tangent.normalized();
    basis.col(2) = bitangent.normalized();
}

template <typename TinyScalar, typename TinyConstants> 
void CollisionBrutal<TinyScalar, TinyConstants>::
LocalCollision() noexcept
{
    // Remove the collision gestion of the parent
    // Build the rotations

    // Contact with an external object
    const unsigned int nbCollisions = m_collisions_infos.size();
    m_local_contact_basis.resize(nbCollisions);

    // Self contact
    const unsigned int nbSelfCollisions = m_self_collisions_infos.size();
    m_local_self_contact_basis.resize(nbSelfCollisions);

    m_remember_self_contact_forces = Eigen::Matrix<TinyScalar, 3, Eigen::Dynamic>::Zero(3, nbSelfCollisions);

// #pragma omp parallel
    {
// #pragma omp for nowait
        for (unsigned int cId = 0u; cId < nbCollisions; ++cId)
        {
            const CollisionInfo<TinyScalar, TinyConstants>& collision_info = m_collisions_infos[cId];
            const Eigen::Matrix<TinyScalar, 3, 1>& normal = collision_info.normal;
            computeBasis(normal, m_local_contact_basis[cId]);
        } // cId

// #pragma omp for
        for (unsigned int scId = 0u; scId < nbSelfCollisions; ++scId)
        {
            const CollisionInfo<TinyScalar, TinyConstants>& self_collision_info = m_self_collisions_infos[scId];
            const Eigen::Matrix<TinyScalar, 3, 1>& normal = self_collision_info.normal;
            computeBasis(normal, m_local_self_contact_basis[scId]);
        } // scId
    }     // omp parallel

//     if (m_handle_self_collision && !m_self_collisions_infos.empty())
//     {
//         TIMER_START(collision_ordering)
//         computeCollisionComputationOrder();
//         const TinyScalar collision_ordering_duration = TIMER_DURATION(collision_ordering, microseconds);
// #ifdef TIMER_PRINT
//         std::cout << "# Collision Ordering: " << collision_ordering_duration << "µs" << std::endl;
// #endif // TIMER_PRINT
//         m_self_collision_ordering_times.push_back(collision_ordering_duration);
//     }
}

template<typename TinyScalar, typename TinyConstants>
Eigen::Matrix<TinyScalar, 3, 1> vec3ToEigen(const Vec3 &x) {
  return Eigen::Matrix<TinyScalar, 3, 1>(
      TinyConstants::scalar_from_double( x[0] ),
      TinyConstants::scalar_from_double( x[1] ),
      TinyConstants::scalar_from_double( x[2] ));
}

template <typename TinyScalar, typename TinyConstants> 
std::vector<CollisionInfo<TinyScalar, TinyConstants>> CollisionBrutal<TinyScalar, TinyConstants>::
getCollisionsInfo(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& positions) const
{
    // TODO: Use better data structures to parallelize.
    std::vector<CollisionInfo<TinyScalar, TinyConstants>> result;

    for (auto impact : obsImpacts) {
      result.push_back(
        CollisionInfo<TinyScalar, TinyConstants>(impact.nodes[0]->index,
            vec3ToEigen<TinyScalar, TinyConstants>(impact.contactPoint),
            vec3ToEigen<TinyScalar, TinyConstants>(impact.n),
            vec3ToEigen<TinyScalar, TinyConstants>(impact.nodes[0]->x - impact.nodes[0]->x0),
            CollisionInfo<TinyScalar, TinyConstants>::default_friction
          ));
    }


    return result;
}


template <typename TinyScalar, typename TinyConstants> 
void CollisionBrutal<TinyScalar, TinyConstants>::
computeSelfCollisionGraph() noexcept
{
    for (std::size_t vertex_index = 0; vertex_index < mCollisionMesh->getNumberOfVertices(); ++vertex_index)
    {
        m_self_collision_graph[vertex_index].clear();
    }

    // Graph Computation
    for (std::size_t self_collision_info_index = 0;
         self_collision_info_index < m_self_collisions_infos.size();
         ++self_collision_info_index)
    {
        const SelfCollisionInfo<TinyScalar, TinyConstants>& self_collision_info =
          m_self_collisions_infos[self_collision_info_index];

        std::size_t target = getSelfCollisionContactPointIndex(self_collision_info);

        std::size_t source = self_collision_info.vertex_index;

        m_self_collision_graph[source].push_back(
          SelfCollisionGraphNeighboor{ target, self_collision_info_index });
        m_self_collision_graph[target].push_back(
          SelfCollisionGraphNeighboor{ source, self_collision_info_index });
    }
}

template <typename TinyScalar, typename TinyConstants> 
std::size_t CollisionBrutal<TinyScalar, TinyConstants>::
getSelfCollisionContactPointIndex(
  const SelfCollisionInfo<TinyScalar, TinyConstants>& self_collision_info) const noexcept
{
    return self_collision_info.vertex_index_b;
}


template <typename TinyScalar, typename TinyConstants> 
void CollisionBrutal<TinyScalar, TinyConstants>::
fillCollisionComputationOrder(const std::vector<std::size_t>& starting_vertices,
                                              std::vector<bool>& visited,
                                              std::size_t initial_computation_index) noexcept
{
    std::vector<std::size_t> visited_last = starting_vertices;
    std::vector<std::size_t> visiting;

    std::size_t computation_index = initial_computation_index;
    while (!visited_last.empty())
    {
        if (m_collision_computation_order.size() <= computation_index)
        {
            m_collision_computation_order.emplace_back();
        }

        for (std::size_t vertex_index : visited_last)
        {
            for (const SelfCollisionGraphNeighboor& neighboor :
                 m_self_collision_graph[vertex_index])
            {
                if (visited[neighboor.collision_index])
                {
                    continue;
                }

                visiting.push_back(neighboor.index);
                visited[neighboor.collision_index] = true;
                m_collision_computation_order[computation_index].push_back(
                  neighboor.collision_index);
            }
        }
        visited_last = std::move(visiting);
        visiting.clear();
        ++computation_index;
    }
}

template <typename TinyScalar, typename TinyConstants> 
bool CollisionBrutal<TinyScalar, TinyConstants>::
vertexIsVisited(std::size_t vertex_index, const std::vector<bool>& visited) const noexcept
{
    for (const auto& neighboor : m_self_collision_graph[vertex_index])
    {
        if (!visited[neighboor.collision_index])
        {
            return false;
        }
    }
    return true;
}


template <typename TinyScalar, typename TinyConstants> 
void CollisionBrutal<TinyScalar, TinyConstants>::
computeCollisionComputationOrder() noexcept
{
    computeSelfCollisionGraph();
    m_collision_computation_order.clear();
    std::vector<bool> visited(m_self_collisions_infos.size(), false);
    std::vector<std::size_t> starting_vertices;

    for (const auto& collision_info : m_collisions_infos)
    {
        starting_vertices.push_back(collision_info.vertex_index);
    }
    fillCollisionComputationOrder(starting_vertices, visited, 0u);

    for (std::size_t vertex_index = 0; vertex_index < mCollisionMesh->getNumberOfVertices(); ++vertex_index)
    {
        if (vertexIsVisited(vertex_index, visited))
        {
            continue;
        }

        if (m_self_collision_graph[vertex_index].size() == 1u)
        {
            starting_vertices.clear();
            starting_vertices.push_back(vertex_index);
            fillCollisionComputationOrder(starting_vertices, visited, 0u);
        }
    }

    for (std::size_t vertex_index = 0; vertex_index < mCollisionMesh->getNumberOfVertices(); ++vertex_index)
    {
        if (vertexIsVisited(vertex_index, visited))
        {
            continue;
        }

        if (m_self_collision_graph[vertex_index].size() > 1u)
        {
            starting_vertices.clear();
            starting_vertices.push_back(vertex_index);
            fillCollisionComputationOrder(starting_vertices, visited, 0u);
        }
    }
}

template <typename TinyScalar, typename TinyConstants> 
std::vector<SelfCollisionInfo<TinyScalar, TinyConstants>> CollisionBrutal<TinyScalar, TinyConstants>::
getSelfCollisionsInfo(TinyScalar tol2) 
{
    std::vector<SelfCollisionInfo<TinyScalar, TinyConstants> > result;
    //filter: each vertex with only one face
    std::map<int, Impact*> closest;
    for (int i = 0; i < impacts.size(); ++i) {
      Impact &imp = impacts[i];
      // std::cout << "impact!";
      // for (int i = 0; i < 4; ++i)
      //   std::cout << " " << imp.nodes[i]->index;
      // std::cout<<std::endl;
      // std::cout << imp.n<<":"<<" ";
      // for (int i = 0; i < 4; ++i)
      //   std::cout << " " << imp.w[i];
      // std::cout << std::endl;
      int vId = imp.nodes[0]->index;
      if (closest.find(vId) == closest.end()) {
        closest[vId] = &imp;
        continue;
      }
      Impact *imp_old = closest[vId];
      if (imp.t < imp_old->t) {
        closest[vId] = &imp;
        continue;
      }
      if (imp_old->t == 1 && imp_old->d > imp.d)
        closest[vId] = &imp;
    }
    //collect: record all vertices with the same face
    for (auto it : closest) {
      int vId = it.first;
      Impact *imp = it.second;
      bool found = false;
      int i;
      for (i = 0; i < result.size(); ++i) {
        bool equal = true;
        for (int k = 0; k < 3; ++k)
          if (result[i].face_indices[k] != imp->nodes[k+1]->index) {
            equal = false;
            break;
          }
        if (equal) {
          found = true;
          break;
        }
      }
      if (found) {
        result[i].vertex_index_list.push_back(vId);
      } else {
        SelfCollisionInfo<TinyScalar, TinyConstants> info;
        info.friction_coefficient = CollisionInfo<TinyScalar, TinyConstants>::default_friction;
        info.normal = vec3ToEigen<TinyScalar, TinyConstants>(imp->n);
        // info.barycentric_coordinates = barycentric_coordinates;            
        info.face_indices = {imp->nodes[1]->index,imp->nodes[2]->index,imp->nodes[3]->index};
        info.vertex_index_list = {vId};
        result.push_back(info);
      }
    }
    // for (int i = 0; i < result.size(); ++i) {
    //   std::cout << "result! (";
    //   for (int k = 0; k < 3; ++k)
    //     std::cout << result[i].face_indices[k] << " ";
    //   std::cout << ") ";
    //   for (int k = 0; k < result[i].vertex_index_list.size(); ++k)
    //     std::cout << result[i].vertex_index_list[k] << " ";
    //   std::cout<<std::endl;
    // }
    return result;

    std::vector<SelfCollisionInfo<TinyScalar, TinyConstants> > collisions_infos;

    std::size_t material_identifier;
    // for (size_t mId = 0u; mId < m_meshes.size(); ++mId)
    // {
    //      const std::shared_ptr<PhysicMesh> mCollisionMesh = m_meshes[mId];
            material_identifier = mCollisionMesh->getMaterialIdentifier();
            collisions_infos = checkMeshSelfCollisions(mCollisionMesh, TinyConstants::sqrt1(tol2));
            for (auto& collision_info : collisions_infos)
            {
                collision_info.friction_coefficient = CollisionInfo<TinyScalar, TinyConstants>::default_friction;
                // collision_info.friction_coefficient = m_friction_coefficients[material_identifier][material_identifier];
            }
            result.insert(result.end(), collisions_infos.begin(), collisions_infos.end());
    // }
    return result;
}






template <typename TinyScalar, typename TinyConstants> 
CGAL::Bbox_3 CollisionBrutal<TinyScalar, TinyConstants>::
getToleranceBoundingBox(const Eigen::Matrix<TinyScalar, 3, 1>& vertex, TinyScalar tolerance)
{

    return CGAL::Bbox_3(TinyConstants::getDouble(vertex.x() - tolerance),
                        TinyConstants::getDouble(vertex.y() - tolerance),
                        TinyConstants::getDouble(vertex.z() - tolerance),
                        TinyConstants::getDouble(vertex.x() + tolerance),
                        TinyConstants::getDouble(vertex.y() + tolerance),
                        TinyConstants::getDouble(vertex.z() + tolerance));
}


template <typename TinyScalar, typename TinyConstants> 
std::vector<std::vector<typename CollisionMesh<TinyScalar, TinyConstants>::Triangle> > 
  CollisionBrutal<TinyScalar, TinyConstants>::
getMeshPotentialTriangleCollision(const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& vertices,
                                  const CollisionMesh<TinyScalar, TinyConstants>* mesh,
                                  TinyScalar tolerance)
{
    std::vector<IndexBox<TinyScalar, TinyConstants> > mesh_boxes;
    for (std::size_t triangle_index = 0; triangle_index < mesh->getNumberOfTriangles(); ++triangle_index)
    {
        IndexBox<TinyScalar, TinyConstants> ib = IndexBox<TinyScalar, TinyConstants>(
          mesh->getTriangleBoundingBox(mesh->getTriangle(triangle_index)), triangle_index);
        mesh_boxes.push_back(ib);
    }

    std::vector<IndexBox<TinyScalar, TinyConstants> > vertices_boxes;
    for (std::size_t vertex_index = 0; vertex_index < vertices.size(); ++vertex_index)
    {
        vertices_boxes.push_back(
          IndexBox<TinyScalar, TinyConstants>(getToleranceBoundingBox(
            vertices[vertex_index], tolerance), vertex_index));
    }

    std::vector<std::vector<typename CollisionMesh<TinyScalar, TinyConstants>::Triangle> > 
    result(vertices.size());

    CGAL::box_intersection_d(
      mesh_boxes.begin(),
      mesh_boxes.end(),
      vertices_boxes.begin(),
      vertices_boxes.end(),
      [&result,&mesh](const IndexBox<TinyScalar, TinyConstants>& triangle_box, 
        const IndexBox<TinyScalar, TinyConstants>& vertex_box) {
          result[vertex_box.m_index].push_back(mesh->getTriangle(triangle_box.m_index));
      });

    return result;
}


template <typename TinyScalar, typename TinyConstants> 
std::vector<std::vector<std::pair<
  Impact*, Eigen::Matrix<TinyScalar, 3, 1>> > > 

CollisionBrutal<TinyScalar, TinyConstants>::
checkMeshAllCollision(const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& vertices,
                      const CollisionMesh<TinyScalar, TinyConstants>* mesh,
                      TinyScalar tolerance)
{
  //   std::vector<std::vector<
    // typename CollisionMesh<TinyScalar, TinyConstants>::Triangle> > collision_to_check =
  //      getMeshPotentialTriangleCollision(vertices, mesh, 1.1 * tolerance);

    std::vector<std::vector<std::pair<
    // typename CollisionMesh<TinyScalar, TinyConstants>::Triangle, Eigen::Matrix<TinyScalar, 3, 1>> > > 
    Impact*, Eigen::Matrix<TinyScalar, 3, 1>> > > 
    result(vertices.size());
    

    for (int i = 0; i < impacts.size(); ++i) {
      Impact &impact = impacts[i];
      std::cout << "impact!";
      for (int i = 0; i < 4; ++i)
        std::cout << " " << impact.nodes[i]->index;
      std::cout<<std::endl;
      std::cout << impact.n<<":"<<" ";
      for (int i = 0; i < 4; ++i)
        std::cout << " " << impact.w[i];
      std::cout << std::endl;
      result[impact.nodes[0]->index].push_back(make_pair(
        &impact,
        // typename CollisionMesh<TinyScalar, TinyConstants>::Triangle{impact.nodes[1]->index,impact.nodes[2]->index,impact.nodes[3]->index},
        Eigen::Matrix<TinyScalar, 3, 1>(-impact.w[1],-impact.w[2],-impact.w[3])
        ));
    }

    // Eigen::Matrix<TinyScalar, 3, 1> current_point;
    // Eigen::Matrix<TinyScalar, 3, 1> current_point_barycentric_coordinate;
    // const TinyScalar tolerance_squared = tolerance * tolerance;
    // for (std::size_t vertex_index = 0; vertex_index < vertices.size(); ++vertex_index)
    // {
    //     for (const typename CollisionMesh<TinyScalar, TinyConstants>::Triangle& 
       //       triangle : collision_to_check[vertex_index])
    //     {
    //         const Eigen::Matrix<TinyScalar, 3, 1>& vertex = vertices[vertex_index];
    //         closestPointToTriangle(vertex,
    //                                mesh->getTriangle(triangle),
    //                                tolerance,
    //                                current_point,
    //                                &current_point_barycentric_coordinate);
            
    //         if ((current_point - vertex).squaredNorm() < tolerance_squared)
    //         {
    //             result[vertex_index].push_back(
    //               std::make_pair(triangle, current_point_barycentric_coordinate));
    //         }
    //     }
    // }
    return result;
}

template <typename TinyScalar, typename TinyConstants> 
std::vector<SelfCollisionInfo<TinyScalar, TinyConstants> >
CollisionBrutal<TinyScalar, TinyConstants>::
checkMeshSelfCollisions(
  const CollisionMesh<TinyScalar, TinyConstants>* mesh, TinyScalar tolerance) noexcept
{
    // std::vector<SelfCollisionInfo<TinyScalar, TinyConstants>> result_me;
    // for (auto impact : impacts) {
    //   std::cout << "impact!";
    //   for (int i = 0; i < 4; ++i)
    //     std::cout << " " << impact.nodes[i]->index;
    //   std::cout<<std::endl;
    //   SelfCollisionInfo<TinyScalar, TinyConstants> self_collision_info;
    //   self_collision_info.normal = vec3ToEigen<TinyScalar, TinyConstants>(impact.n);
    //   Eigen::Matrix<TinyScalar, 3, 1> barycentric_coordinates;
    //   for (int k = 0; k < 3; ++k)
    //     barycentric_coordinates[k] = -impact.w[k+1];
    //   self_collision_info.barycentric_coordinates = barycentric_coordinates;            
    //   self_collision_info.face_indices = {impact.nodes[1]->index,impact.nodes[2]->index,impact.nodes[3]->index};
    //   self_collision_info.vertex_index_b = impact.nodes[1]->index;
    //   self_collision_info.contact_point = vec3ToEigen<TinyScalar, TinyConstants>(impact.contactPoint);
    //   self_collision_info.vertex_index = impact.nodes[0]->index;
    //   result_me.push_back(self_collision_info);
    // }
    // return result_me;

    std::vector<Eigen::Matrix<TinyScalar, 3, 1>> vertices = mesh->getVertices();

  

    std::vector<std::vector<std::pair<
    // typename CollisionMesh<TinyScalar, TinyConstants>::Triangle, Eigen::Matrix<TinyScalar, 3, 1>> > > 
    Impact*, Eigen::Matrix<TinyScalar, 3, 1>> > > 
    collisions_infos = checkMeshAllCollision(vertices, mesh, tolerance);


    std::vector<SelfCollisionInfo<TinyScalar, TinyConstants>> result;
    typename CollisionMesh<TinyScalar, TinyConstants>::Triangle triangle;
    Eigen::Matrix<TinyScalar, 3, 1> barycentric_coordinates;
    Eigen::Matrix<TinyScalar, 3, 1> contact_point;
    SelfCollisionInfo<TinyScalar, TinyConstants> self_collision_info;
    std::array<bool, 3> face_vertices_candidate;

    std::vector<std::vector<size_t> > registered;
    registered.resize(vertices.size());
    Impact *impact;
    
    // Since we only do node-node contact, we pair each vertex in collision with
    // one other vertex. For each face, we make sure that the vertex is paired
    // with only one vertex of this face.  We will call the vertex of the
    // outside loop the colliding vertex.
    for (std::size_t vertex_index = 0; vertex_index < vertices.size(); ++vertex_index)
    {
        for (const auto& collision_info : collisions_infos[vertex_index])
        {
            boost::tie(impact, barycentric_coordinates) = collision_info;
            triangle = 
        typename CollisionMesh<TinyScalar, TinyConstants>::Triangle{impact->nodes[1]->index,impact->nodes[2]->index,impact->nodes[3]->index};
            // If the colliding vertex is part of the triangle or if it is on
            // the edge we don't count this as a collision.
            // if (std::find(triangle.vertex_indices.begin(),
            //               triangle.vertex_indices.end(),
            //               vertex_index) != triangle.vertex_indices.end() ||
            //     barycentric_coordinates.x() == 0 || barycentric_coordinates.y() == 0 ||
            //     barycentric_coordinates.z() == 0)
            // {

            //     continue;
            // }

            // The colliding vertex is already paired with a vertex of this
            // triangle so we cannot count a collision with this triangle.
            // if (isTriangleVertexIndexIn(triangle, registered[vertex_index]) < 3u)
            // {
              
            //     continue;
            // }

            // We try to find a vertex within the triangle that is not is
            // contact with a vertex adjacent to the colliding vertex. If we
            // were to link a vertex not satisfying this conditiong with the
            // colliding vertex. The linked vertex would be link to two vertex
            // that are part of the same triangle, the colliding vertex and its
            // adjacent vertex.
            std::fill(face_vertices_candidate.begin(), face_vertices_candidate.end(), true);
            // for (std::size_t one_ring_vertex_index : mesh->getVerticesAroundVertexRange(vertex_index))
            // {
            //     std::size_t triangle_vertex_index = isTriangleVertexIndexIn(triangle, registered[one_ring_vertex_index]);
            //     if (triangle_vertex_index < 3u)
            //     {
            //         face_vertices_candidate[triangle_vertex_index] = false;
            //     }
            // }

            // Among the candidate obtained in the previous loop. We take the
            // closest one.
            std::size_t closest_candidate_index = 3u;
            TinyScalar closest_candidate_barycentric_component = 0;
            for (std::size_t candidate_index = 0; candidate_index < 3u; ++candidate_index)
            {
                if (!face_vertices_candidate[candidate_index])
                {
            
                    continue;
                }

                if (barycentric_coordinates[candidate_index] > closest_candidate_barycentric_component)
                {
                    closest_candidate_index = candidate_index;
                    closest_candidate_barycentric_component = barycentric_coordinates[candidate_index];
                }
            }

            //All the triangle vertices were already paired with one vertex of the one ring. Therefore, we don't pair this one.
            if (closest_candidate_index >= 3u)
            {
                
                continue;
            }

            contact_point = mesh->getTriangle(triangle) * barycentric_coordinates;
            self_collision_info.normal = vec3ToEigen<TinyScalar, TinyConstants>(impact->n);
              // (mesh->getVertex(vertex_index) - contact_point).normalized();
            
            std::size_t paired_vertex_index = triangle.vertex_indices[closest_candidate_index];
            // barycentric_coordinates[closest_candidate_index] = 1.;
            // barycentric_coordinates[(closest_candidate_index + 1) % 3] = 0.;
            // barycentric_coordinates[(closest_candidate_index + 2) % 3] = 0.;

            // Keep track of which vertex has been paired with which.
            registered[paired_vertex_index].push_back(vertex_index);
            registered[vertex_index].push_back(paired_vertex_index);

            contact_point = mesh->getTriangle(triangle) * barycentric_coordinates;
            self_collision_info.barycentric_coordinates = barycentric_coordinates;            
            self_collision_info.face_indices = triangle.vertex_indices;
            self_collision_info.vertex_index_b = triangle.vertex_indices[closest_candidate_index];
            self_collision_info.contact_point = contact_point;
// #warning Suppose there is only one mesh simulated
            self_collision_info.vertex_index = vertex_index;
            result.push_back(std::move(self_collision_info));
              
            }
            
            
        }

    return result;
}


/**
 * @brief From Christer Ericson -- Real-Time Collision Detection (p141)
 *        (From an implementation by Gilles Daviet)
 *
 * @param p     Vertex to project
 * @param triangle Matrix whose column are the positions of the triangle vertices
 * @param tol2  Tolerance for the collision
 *
 * @return The barycentric coordinate of the closest point to p inside (abc)
 */

template <typename TinyScalar, typename TinyConstants> 
void CollisionBrutal<TinyScalar, TinyConstants>::
closestPointToTriangle(const Eigen::Matrix<TinyScalar, 3, 1>& p,
                       const Eigen::Matrix<TinyScalar, 3, 3>& triangle,
                       TinyScalar tol2,
                       Eigen::Matrix<TinyScalar, 3, 1>& closest_point,
                       Eigen::Matrix<TinyScalar, 3, 1>* barycentric_coordinates_ptr) noexcept
{
    const auto a = triangle.col(0);
    const auto b = triangle.col(1);
    const auto c = triangle.col(2);

    // Check if P in vertex region outside A
    const Eigen::Matrix<TinyScalar, 3, 1> ap = p - a;
    const Eigen::Matrix<TinyScalar, 3, 1> ab = b - a;
    const Eigen::Matrix<TinyScalar, 3, 1> ac = c - a;

    // Early exit
    // Todo : tune tolerance

    if (ap.squaredNorm() - ab.squaredNorm() - ac.squaredNorm() > tol2)
    {
        closest_point = a;
        if (barycentric_coordinates_ptr)
        {
            *barycentric_coordinates_ptr = Eigen::Matrix<TinyScalar, 3, 1>{ 1, 0, 0 };
        }
        return;
    }
    const TinyScalar d1 = ab.dot(ap);
    const TinyScalar d2 = ac.dot(ap);
    if (d1 <= 0.0d && d2 <= 0.0d)
    {
        closest_point = a;
        if (barycentric_coordinates_ptr)
        {
            *barycentric_coordinates_ptr = Eigen::Matrix<TinyScalar, 3, 1>{ 1, 0, 0 };
        }
        return;
    }
    // Check if P in vertex region outside B
    const Eigen::Matrix<TinyScalar, 3, 1> bp = p - b;
    const TinyScalar d3 = ab.dot(bp);
    const TinyScalar d4 = ac.dot(bp);
    if (d3 >= 0.0d && d4 <= d3)
    {
        closest_point = b;
        if (barycentric_coordinates_ptr)
        {
            *barycentric_coordinates_ptr = Eigen::Matrix<TinyScalar, 3, 1>{ 0, 1, 0 };
        }
        return;
    }
    // Check if P in edge region of AB, if so return projection of P onto AB
    const TinyScalar vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0d && d1 >= 0.0d && d3 <= 0.0d)
    {
        const TinyScalar v = d1 / (d1 - d3);
        closest_point = a + v * ab;
        if (barycentric_coordinates_ptr)
        {
            *barycentric_coordinates_ptr = Eigen::Matrix<TinyScalar, 3, 1>{ 1 - v, v, 0 };
        }
        return;
    }
    // Check if P in vertex region outside C
    const Eigen::Matrix<TinyScalar, 3, 1> cp = p - c;
    const TinyScalar d5 = ab.dot(cp);
    const TinyScalar d6 = ac.dot(cp);
    if (d6 >= 0.0d && d5 <= d6)
    {
        closest_point = c;
        if (barycentric_coordinates_ptr)
        {
            *barycentric_coordinates_ptr = Eigen::Matrix<TinyScalar, 3, 1>{ 0, 0, 1 };
        }
        return;
    }
    // Check if P in edge region of AC, if so return projection of P onto AC
    const TinyScalar vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0d && d2 >= 0.0d && d6 <= 0.0d)
    {
        const TinyScalar w = d2 / (d2 - d6);
        closest_point = a + w * ac;
        if (barycentric_coordinates_ptr)
        {
            *barycentric_coordinates_ptr = Eigen::Matrix<TinyScalar, 3, 1>{ 1 - w, 0, w };
        }
        return;
    }
    // Check if P in edge region of BC, if so return projection of P onto BC
    TinyScalar va = d3 * d6 - d5 * d4;
    if (va <= 0.0d && (d4 - d3) >= 0.0d && (d5 - d6) >= 0.0d)
    {
        const TinyScalar w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        closest_point = b + w * (c - b);
        if (barycentric_coordinates_ptr)
        {
            *barycentric_coordinates_ptr = Eigen::Matrix<TinyScalar, 3, 1>{ 0, 1 - w, w };
        }
        return;
    }
    // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
    TinyScalar denom = 1.0d / (va + vb + vc);
    TinyScalar v = vb * denom;
    TinyScalar w = vc * denom;
    closest_point = a + ab * v + ac * w;
    if (barycentric_coordinates_ptr)
    {
        *barycentric_coordinates_ptr = Eigen::Matrix<TinyScalar, 3, 1>{ 1. - v - w, v, w };
    }
}


/**
 *  @return i where triangle.vertices_indices[i] is in vertices_indices. If there is no such i, returns triange.vertices_indices.size().
 */
template <typename TinyScalar, typename TinyConstants> 
std::size_t CollisionBrutal<TinyScalar, TinyConstants>::
isTriangleVertexIndexIn(
  const typename CollisionMesh<TinyScalar, TinyConstants>::Triangle& triangle, 
  const std::vector<size_t>& vertices_indices)
{
  auto triangle_vertex_index_it = std::find_first_of(
    triangle.vertex_indices.begin(),
    triangle.vertex_indices.end(),
    vertices_indices.begin(),
    vertices_indices.end());
  return std::distance(triangle.vertex_indices.begin(), triangle_vertex_index_it);
}


}
#endif