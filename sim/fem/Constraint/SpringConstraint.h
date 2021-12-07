#ifndef __SPRING_CONSTRAINT_H__
#define __SPRING_CONSTRAINT_H__ 
#include "Constraint.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>

namespace FEM
{
enum ConstraintType;

template <typename TinyScalar, typename TinyConstants> 
class SpringConstraint : public Constraint<TinyScalar, TinyConstants>
{
public:
  SpringConstraint(int i0, int i1, TinyScalar rest_length, TinyScalar weight);

  int GetDof() override;
  ConstraintType GetType() override;

  void  EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets);
  void  EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets);
  void  EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
  void  GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d);

  int&       GetI0() {return mi0;}
  Eigen::Matrix<TinyScalar, 3, 1>& GetP()  {return mp;}
  TinyScalar mStretch;
	int 			mi0,mi1;
  TinyScalar    m_rest_length;

protected:
  Eigen::Matrix<TinyScalar, 3, 1> mp;
  Eigen::Matrix<TinyScalar, 3, 1> md;
};




template <typename TinyScalar, typename TinyConstants> 
SpringConstraint<TinyScalar, TinyConstants>::
SpringConstraint(int i0, int i1, TinyScalar rest_length, TinyScalar weight)
  :Constraint<TinyScalar, TinyConstants>(weight), mi0(i0), 
  mi1(i1), m_rest_length(rest_length), mStretch(sqrt(weight))
{
  md.setZero();
}

template <typename TinyScalar, typename TinyConstants> 
int SpringConstraint<TinyScalar, TinyConstants>::
GetDof()
{
  return 1;
}

template <typename TinyScalar, typename TinyConstants> 
ConstraintType SpringConstraint<TinyScalar, TinyConstants>::
GetType()    
{
  return ConstraintType::ATTACHMENT; 
}

template <typename TinyScalar, typename TinyConstants> 
void SpringConstraint<TinyScalar, TinyConstants>::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets)
{
  Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> Ai(1*3,3*2);

	Ai<<
		1,0,0,-1,0,0,
    0,1,0,0,-1,0,
    0,0,1,0,0,-1;
    
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> MuAiT = mStretch*Ai.transpose();
	int idx[4] = {mi0,mi1};
	for(int i =0;i<2;i++)
	{
    int j = 0;
    J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+0, 3*(index+j)+0, MuAiT(3*i+0, 3*j+0)));
    J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+0, 3*(index+j)+1, MuAiT(3*i+0, 3*j+1)));
    J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+0, 3*(index+j)+2, MuAiT(3*i+0, 3*j+2)));
    J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+1, 3*(index+j)+0, MuAiT(3*i+1, 3*j+0)));
    J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+1, 3*(index+j)+1, MuAiT(3*i+1, 3*j+1)));
    J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+1, 3*(index+j)+2, MuAiT(3*i+1, 3*j+2)));
    J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+2, 3*(index+j)+0, MuAiT(3*i+2, 3*j+0)));
    J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+2, 3*(index+j)+1, MuAiT(3*i+2, 3*j+1)));
    J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+2, 3*(index+j)+2, MuAiT(3*i+2, 3*j+2)));
	}
}

template <typename TinyScalar, typename TinyConstants> 
void SpringConstraint<TinyScalar, TinyConstants>::
EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets)
{
  Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> Ai(1*3,3*2);

	Ai<<
		1,0,0,-1,0,0,
    0,1,0,0,-1,0,
    0,0,1,0,0,-1;
    
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> MuAiTAi = mStretch*mStretch*((Ai.transpose())*Ai);
	int idx[2] = {mi0,mi1};
	//MuAiT --- 12x12 matrix
	for(int i =0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			//MuAiTAi.block [i,j] -- 3x3 matrix
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+0, 3*idx[j]+0, MuAiTAi(3*i+0, 3*j+0)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+0, 3*idx[j]+1, MuAiTAi(3*i+0, 3*j+1)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+0, 3*idx[j]+2, MuAiTAi(3*i+0, 3*j+2)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+1, 3*idx[j]+0, MuAiTAi(3*i+1, 3*j+0)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+1, 3*idx[j]+1, MuAiTAi(3*i+1, 3*j+1)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+1, 3*idx[j]+2, MuAiTAi(3*i+1, 3*j+2)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+2, 3*idx[j]+0, MuAiTAi(3*i+2, 3*j+0)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+2, 3*idx[j]+1, MuAiTAi(3*i+2, 3*j+1)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+2, 3*idx[j]+2, MuAiTAi(3*i+2, 3*j+2)));
		}
	}
}

template <typename TinyScalar, typename TinyConstants> 
void SpringConstraint<TinyScalar, TinyConstants>::
EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x)
{
  Eigen::Matrix<TinyScalar, 3, 1> x0(x.template block<3,1>(mi0*3,0));
  Eigen::Matrix<TinyScalar, 3, 1> x1(x.template block<3,1>(mi1*3,0));
  Eigen::Matrix<TinyScalar, 3, 1> e = x0 - x1;
  
  // std::cout << x0 << " p1\n";
  // std::cout << x1 << " p2\n";
  // std::cout << mStretch << " m_weight\n";
  // std::cout << m_rest_length << " m_rest_length\n";
  e.normalize();  
  md = mStretch * m_rest_length * e;
}

template <typename TinyScalar, typename TinyConstants> 
void SpringConstraint<TinyScalar, TinyConstants>::
GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d)
{
  d.template block<3,1>(index*3,0) = md;
  index++;
}
};
#endif