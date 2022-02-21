#ifndef __ATTACHMENT_CONSTRAINT_H__
#define __ATTACHMENT_CONSTRAINT_H__	
#include "Constraint.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>

namespace FEM
{
enum ConstraintType;

template <typename TinyScalar, typename TinyConstants> 
class AttachmentConstraint : public Constraint<TinyScalar, TinyConstants>
{
public:
	AttachmentConstraint(const TinyScalar& stiffness,int i0,const Eigen::Matrix<TinyScalar, 3, 1>& p);

	int GetDof() override;
	ConstraintType GetType() override;

	void	EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets);
	void	EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets);
	void 	EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
	void 	GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d);
	void	fixIndex(int offset);

	int&			 GetI0() {return mi0;}
	Eigen::Matrix<TinyScalar, 3, 1>& GetP()  {return mp;}

	TinyScalar mStiffness;
// protected:
	int mi0;
	Eigen::Matrix<TinyScalar, 3, 1> mp;
	Eigen::Matrix<TinyScalar, 3, 1> md;
};




template <typename TinyScalar, typename TinyConstants> 
AttachmentConstraint<TinyScalar, TinyConstants>::
AttachmentConstraint(const TinyScalar& stiffness,int i0,const Eigen::Matrix<TinyScalar, 3, 1>& p)
	:Constraint<TinyScalar, TinyConstants>(stiffness),mi0(i0),mp(p)
{
	mStiffness = stiffness;
	md.setZero();
}

template <typename TinyScalar, typename TinyConstants> 
int AttachmentConstraint<TinyScalar, TinyConstants>::
GetDof()
{
	return 1;
}

template <typename TinyScalar, typename TinyConstants> 
ConstraintType AttachmentConstraint<TinyScalar, TinyConstants>::
GetType()	   
{
	return ConstraintType::ATTACHMENT; 
}

template <typename TinyScalar, typename TinyConstants> 
void AttachmentConstraint<TinyScalar, TinyConstants>::fixIndex(int offset) {
  mi0 += offset;
}

template <typename TinyScalar, typename TinyConstants> 
void AttachmentConstraint<TinyScalar, TinyConstants>::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets)
{
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+0,3*index+0,mStiffness));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+1,3*index+1,mStiffness));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+2,3*index+2,mStiffness));
}

template <typename TinyScalar, typename TinyConstants> 
void AttachmentConstraint<TinyScalar, TinyConstants>::
EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets)
{
	L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+0,3*mi0+0,mStiffness));
	L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+1,3*mi0+1,mStiffness));
	L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+2,3*mi0+2,mStiffness));
}

template <typename TinyScalar, typename TinyConstants> 
void AttachmentConstraint<TinyScalar, TinyConstants>::
EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x)
{
	md = mp;
}

template <typename TinyScalar, typename TinyConstants> 
void AttachmentConstraint<TinyScalar, TinyConstants>::
GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d)
{
	d.template block<3,1>(index*3,0) = md;
	index++;
}
};
#endif