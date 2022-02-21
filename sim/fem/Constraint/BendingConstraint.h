#ifndef __BENDING_CONSTRAINT_H__
#define __BENDING_CONSTRAINT_H__ 
#include "Constraint.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>

namespace FEM
{
enum ConstraintType;

template <typename TinyScalar, typename TinyConstants> 
class BendingConstraint : public Constraint<TinyScalar, TinyConstants>
{
public:
  BendingConstraint(const std::vector<int>& vertex_indices,
                  const Eigen::Matrix<TinyScalar, 3, 1>& middle_vertex_rest_state,
                  const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& rest_state,
                  TinyScalar weight);

  int GetDof() override;
  ConstraintType GetType() override;

  void  EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets);
  void  EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets);
  void  EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
  void  GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d);

  void  fixIndex(int offset);
  int&       GetI0() {return mi0;}
  Eigen::Matrix<TinyScalar, 3, 1>& GetP()  {return mp;}

  TinyScalar mBend;
  std::vector<int> m_relevant_indices;  
  /// @brief Norm of the Laplacian of the rest shape
  TinyScalar m_rest_state_laplace_beltrami_norm;
  /// @brief Normalization coefficient for the Laplacian-Beltrami
  TinyScalar m_laplace_beltrami_coefficient_sum;
  std::vector<TinyScalar> m_beltrami_coefficients;
  /// @brief A in the PD potentiels
  Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> m_A;
  Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> m_ATA;


  Eigen::Matrix<TinyScalar, 3, 1> applyA(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& relevant_vector) const noexcept;
  Eigen::Matrix<TinyScalar, 3, Eigen::Dynamic> getLaplaceBeltramiDicretisation(const Eigen::Matrix<TinyScalar, 3, 1>& middle_vertex,
                                  const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& ring_vertices);
  std::vector<TinyScalar> getLaplaceBeltramiDiscretisationMeanValueCoefficients(
    const Eigen::Matrix<TinyScalar, 3, 1>& middle_vertex,
    const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& ring_vertices);
  TinyScalar getMeanValue(const Eigen::Matrix<TinyScalar, 3, 1>& v1, const Eigen::Matrix<TinyScalar, 3, 1>& v2, const Eigen::Matrix<TinyScalar, 3, 1>& v3);
  TinyScalar getSine(const Eigen::Matrix<TinyScalar, 3, 1>& v1, const Eigen::Matrix<TinyScalar, 3, 1>& v2, const Eigen::Matrix<TinyScalar, 3, 1>& v3);
  TinyScalar getCosine(const Eigen::Matrix<TinyScalar, 3, 1>& v1, const Eigen::Matrix<TinyScalar, 3, 1>& v2, const Eigen::Matrix<TinyScalar, 3, 1>& v3);

protected:
  int mi0;
  Eigen::Matrix<TinyScalar, 3, 1> mp;
  Eigen::Matrix<TinyScalar, 3, 1> md;
  std::vector<int> mRelevantIndices;
};





template <typename TinyScalar, typename TinyConstants> 
BendingConstraint<TinyScalar, TinyConstants>::
BendingConstraint(const std::vector<int>& vertex_indices,
                  const Eigen::Matrix<TinyScalar, 3, 1>& middle_vertex_rest_state,
                  const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& rest_state,
                  TinyScalar weight)
  :Constraint<TinyScalar, TinyConstants>(weight),
  mBend(weight), m_relevant_indices(vertex_indices)
{
    assert(rest_state.size() > 3);
    assert(vertex_indices.size() == rest_state.size() + 1);
    for (const auto& idx : vertex_indices)
      mRelevantIndices.push_back(idx);
    // TODO: factorize with update()
    m_A = weight * getLaplaceBeltramiDicretisation(middle_vertex_rest_state, rest_state);
    m_ATA = m_A.transpose() * m_A;

    // This equality come from the definition of the discretisation of
    // the laplace beltrami operator.
    m_laplace_beltrami_coefficient_sum = -m_A(0, 0);

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> generalized_rest_state((rest_state.size() + 1) * 3);
    for (std::size_t i = 0; i < rest_state.size(); ++i)
    {
      generalized_rest_state.segment<3>(3 * (i + 1)) = rest_state[i];
    }
    generalized_rest_state.segment<3>(3 * 0) = middle_vertex_rest_state;

    m_rest_state_laplace_beltrami_norm = applyA(generalized_rest_state).norm() ; //(m_A * generalized_rest_state).norm();

    const std::vector<TinyScalar> coeffs =
      getLaplaceBeltramiDiscretisationMeanValueCoefficients(middle_vertex_rest_state, rest_state);
    m_beltrami_coefficients.reserve(m_relevant_indices.size());
    m_beltrami_coefficients.push_back(0.);
    for (size_t i = 0u; i < coeffs.size(); ++i)
    {
      m_beltrami_coefficients[0] += coeffs[i];
      m_beltrami_coefficients.push_back(-coeffs[i]);
    }
}

template <typename TinyScalar, typename TinyConstants> 
int BendingConstraint<TinyScalar, TinyConstants>::
GetDof()
{
  return 1;
}

template <typename TinyScalar, typename TinyConstants> 
ConstraintType BendingConstraint<TinyScalar, TinyConstants>::
GetType()    
{
  return ConstraintType::ATTACHMENT; 
}

template <typename TinyScalar, typename TinyConstants> 
void BendingConstraint<TinyScalar, TinyConstants>::fixIndex(int offset) {
  mi0 += offset;
}
template <typename TinyScalar, typename TinyConstants> 
void BendingConstraint<TinyScalar, TinyConstants>::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets)
{

  for(int i = 0;i < mRelevantIndices.size(); i++) {
    J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mRelevantIndices[i]+0, 3*index+0, mBend *  m_beltrami_coefficients[i]));
    J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mRelevantIndices[i]+1, 3*index+1, mBend *  m_beltrami_coefficients[i]));
    J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mRelevantIndices[i]+2, 3*index+2, mBend *  m_beltrami_coefficients[i]));
  }
}

template <typename TinyScalar, typename TinyConstants> 
void BendingConstraint<TinyScalar, TinyConstants>::
EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets)
{

  Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> MuAiTAi = mBend*((m_A.transpose())*m_A);

  for (size_t i = 0u; i < m_relevant_indices.size(); ++i)
  {
      for (size_t j = 0u; j < m_relevant_indices.size(); ++j)
      {
          //MuAiTAi.block [i,j] -- 3x3 matrix
          L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mRelevantIndices[i]+0, 3*mRelevantIndices[j]+0, MuAiTAi(3*i+0, 3*j+0)));
          L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mRelevantIndices[i]+0, 3*mRelevantIndices[j]+1, MuAiTAi(3*i+0, 3*j+1)));
          L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mRelevantIndices[i]+0, 3*mRelevantIndices[j]+2, MuAiTAi(3*i+0, 3*j+2)));
          L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mRelevantIndices[i]+1, 3*mRelevantIndices[j]+0, MuAiTAi(3*i+1, 3*j+0)));
          L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mRelevantIndices[i]+1, 3*mRelevantIndices[j]+1, MuAiTAi(3*i+1, 3*j+1)));
          L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mRelevantIndices[i]+1, 3*mRelevantIndices[j]+2, MuAiTAi(3*i+1, 3*j+2)));
          L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mRelevantIndices[i]+2, 3*mRelevantIndices[j]+0, MuAiTAi(3*i+2, 3*j+0)));
          L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mRelevantIndices[i]+2, 3*mRelevantIndices[j]+1, MuAiTAi(3*i+2, 3*j+1)));
          L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mRelevantIndices[i]+2, 3*mRelevantIndices[j]+2, MuAiTAi(3*i+2, 3*j+2)));
      }
  }
}

template <typename TinyScalar, typename TinyConstants> 
void BendingConstraint<TinyScalar, TinyConstants>::
EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x)
{    
  Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> relevant_vector(m_relevant_indices.size() * 3);
  for (std::size_t i = 0; i < m_relevant_indices.size(); ++i)
  {
    relevant_vector.segment<3>(3 * i) =
      x.segment<3>(3 * m_relevant_indices[i]);
  }
  
  Eigen::Matrix<TinyScalar, 3, 1> laplace_beltrami = applyA(relevant_vector);

  TinyScalar laplace_beltrami_norm = laplace_beltrami.norm();

  // If the laplace beltrami is the null vector we do not normalize it
  if (laplace_beltrami_norm > 1e-25)
  {
    laplace_beltrami *= m_rest_state_laplace_beltrami_norm / laplace_beltrami_norm;
    md = mBend * laplace_beltrami;
  }
  else
  {
    md.setZero();
  }
}

template <typename TinyScalar, typename TinyConstants> 
void BendingConstraint<TinyScalar, TinyConstants>::
GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d)
{
  d.template block<3,1>(index*3,0) = md;
  index++;
}




template <typename TinyScalar, typename TinyConstants> 
TinyScalar BendingConstraint<TinyScalar, TinyConstants>::
getMeanValue(const Eigen::Matrix<TinyScalar, 3, 1>& v1, const Eigen::Matrix<TinyScalar, 3, 1>& v2, const Eigen::Matrix<TinyScalar, 3, 1>& v3)
{
    TinyScalar numerator = (1 - getCosine(v1, v2, v3));
    if (numerator == 0)
    {
        return 0;
    }
    return numerator / getSine(v1, v2, v3);
}

template <typename TinyScalar, typename TinyConstants> 
std::vector<TinyScalar> BendingConstraint<TinyScalar, TinyConstants>::
getLaplaceBeltramiDiscretisationMeanValueCoefficients(
  const Eigen::Matrix<TinyScalar, 3, 1>& middle_vertex,
  const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& ring_vertices)
{
    // There can't be only two triangle around the middle vertex.
    // It would mean that one of the triangle is flat, and therefore the
    // triangulation of the mesh would be wrong. Or the vertex would be on the
    // border and therefore the Laplace-Beltrami does not exist for this point.
    assert(ring_vertices.size() >= 3);

    TinyScalar clockwise_mean_value, counter_clockwise_mean_value;
    TinyScalar distance;
    std::vector<TinyScalar> result;
    for (std::size_t i = 0; i < ring_vertices.size(); ++i)
    {
        const Eigen::Matrix<TinyScalar, 3, 1>& vertex = ring_vertices[i];
        const Eigen::Matrix<TinyScalar, 3, 1>& clockwise_vertex = ring_vertices[(i + 1) % ring_vertices.size()];
        const Eigen::Matrix<TinyScalar, 3, 1>& counter_clockwise_vertex =
          ring_vertices[(i + ring_vertices.size() - 1) % ring_vertices.size()];

        // We can't have a flat triangle around the middle vertex. If that
        // was the case, the triangulation of the mesh would be wrong or
        // the middle vertex would be on a border.
        assert(!areCollinear(vertex, clockwise_vertex, middle_vertex));
        assert(!areCollinear(vertex, middle_vertex, clockwise_vertex));

        distance = (vertex - middle_vertex).norm();

        clockwise_mean_value = getMeanValue(vertex, middle_vertex, clockwise_vertex);

        counter_clockwise_mean_value =
          getMeanValue(vertex, counter_clockwise_vertex, middle_vertex);

        result.push_back((clockwise_mean_value + counter_clockwise_mean_value) / distance);
    }
    return result;
}

template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, 3, Eigen::Dynamic> BendingConstraint<TinyScalar, TinyConstants>::
getLaplaceBeltramiDicretisation(const Eigen::Matrix<TinyScalar, 3, 1>& middle_vertex,
                                const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& ring_vertices)
{
    // There can't be only two triangle around the middle vertex.
    // It would mean that one of the triangle is flat, and therefore the
    // triangulation of the mesh would be wrong. Or the vertex would be on the
    // border and therefore the Laplace-Beltrami does not exist for this point.
    assert(ring_vertices.size() >= 3);

    Eigen::Matrix<TinyScalar, 3, Eigen::Dynamic> matrix = Eigen::Matrix<TinyScalar, 3, Eigen::Dynamic>::Zero(3, (ring_vertices.size() + 1) * 3);

    std::vector<TinyScalar> coeff =
      getLaplaceBeltramiDiscretisationMeanValueCoefficients(middle_vertex, ring_vertices);

    Eigen::Matrix<TinyScalar, 3, 3> middle_vertex_block = Eigen::Matrix<TinyScalar, 3, 3>::Zero(3, 3);
    for (std::size_t i = 0; i < ring_vertices.size(); ++i)
    {
        matrix.template block<3, 3>(0, (1 + i) * 3) = -coeff[i] * Eigen::Matrix<TinyScalar, 3, 3>::Identity();
        middle_vertex_block += coeff[i] * Eigen::Matrix<TinyScalar, 3, 3>::Identity();
    }
    matrix.template block<3, 3>(0, 0) = middle_vertex_block;

    return matrix;
}

template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, 3, 1> BendingConstraint<TinyScalar, TinyConstants>::
applyA(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& relevant_vector) const noexcept
{
    Eigen::Matrix<TinyScalar, 3, 1> result(Eigen::Matrix<TinyScalar, 3, 1>::Zero());
    for (std::size_t vertex_index = 0; vertex_index < m_beltrami_coefficients.size();
         ++vertex_index)
    {
        result += 
          m_beltrami_coefficients[vertex_index] * relevant_vector.segment<3>(3 * vertex_index);
    }
    return result;
}


template <typename TinyScalar, typename TinyConstants> 
TinyScalar BendingConstraint<TinyScalar, TinyConstants>::
getSine(const Eigen::Matrix<TinyScalar, 3, 1>& v1, const Eigen::Matrix<TinyScalar, 3, 1>& v2, const Eigen::Matrix<TinyScalar, 3, 1>& v3)
{
    Eigen::Matrix<TinyScalar, 3, 1> edge1 = v1 - v2;
    Eigen::Matrix<TinyScalar, 3, 1> edge2 = v3 - v2;
    return edge1.cross(edge2).norm() / (edge1.norm() * edge2.norm());
}

template <typename TinyScalar, typename TinyConstants> 
TinyScalar BendingConstraint<TinyScalar, TinyConstants>::
getCosine(const Eigen::Matrix<TinyScalar, 3, 1>& v1, const Eigen::Matrix<TinyScalar, 3, 1>& v2, const Eigen::Matrix<TinyScalar, 3, 1>& v3)
{
    Eigen::Matrix<TinyScalar, 3, 1> edge1 = v1 - v2;
    Eigen::Matrix<TinyScalar, 3, 1> edge2 = v3 - v2;
    return edge1.dot(edge2) / (edge1.norm() * edge2.norm());
}


};
#endif