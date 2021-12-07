#ifndef __TRIANGLE_MUSCLE_CONSTRAINT_H__
#define __TRIANGLE_MUSCLE_CONSTRAINT_H__
#include "Constraint.h"
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/Sparse>
#include <iostream>
namespace Eigen
{
using Matrix32d = Matrix<TinyScalar, 3, 2>;
};
namespace FEM
{
template <typename TinyScalar, typename TinyConstants> 
class TriangleMuscleConstraint : public Constraint<TinyScalar, TinyConstants>
{
public:
	TriangleMuscleConstraint(TinyScalar stiffness,const Eigen::Matrix<TinyScalar, 2, 1>& fiber_direction,int i0,int i1,int i2,TinyScalar area,const Eigen::Matrix2d& invDm);
	int GetDof() override;

	void EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets);
	void EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets);
	void EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
	void GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d);

	int GetI0() {return mi0;}
	int GetI1() {return mi1;}
	int GetI2() {return mi2;}

	void SetActivationLevel(TinyScalar a){mActivationLevel = std::min(1.0,std::max(0.0,a));}
	TinyScalar GetActivationLevel(){return mActivationLevel;}
	const Eigen::Matrix<TinyScalar, 2, 1>& GetFiberDirection() {return mFiberDirection;}
protected:
	Eigen::Matrix<TinyScalar, 2, 1> mFiberDirection;
	TinyScalar			mActivationLevel;

	int mi0,mi1,mi2;
	TinyScalar mArea;
	Eigen::Matrix2d mInvDm;

	Eigen::Matrix<TinyScalar, 3, 1> md;
};




template <typename TinyScalar, typename TinyConstants> 
TriangleMuscleConstraint<TinyScalar, TinyConstants>::
TriangleMuscleConstraint(TinyScalar stiffness,const Eigen::Matrix<TinyScalar, 2, 1>& fiber_direction,int i0,int i1,int i2,TinyScalar area,const Eigen::Matrix2d& invDm)
	:Constraint<TinyScalar, TinyConstants>(stiffness),mFiberDirection(fiber_direction),
	mi0(i0),mi1(i1),mi2(i2),mArea(area),mInvDm(invDm),mActivationLevel(0.0)
{

}

template <typename TinyScalar, typename TinyConstants> 
int TriangleMuscleConstraint<TinyScalar, TinyConstants>::
GetDof()
{
	return 1;
}


template <typename TinyScalar, typename TinyConstants> 
void TriangleMuscleConstraint<TinyScalar, TinyConstants>::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets)
{
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> Ai(3,9);

	Eigen::Matrix<TinyScalar, 2, 1> v = mInvDm*mFiberDirection;

	Ai<<
		-v[0]-v[1],0,0,v[0],0,0,v[1],0,0,
		0,-v[0]-v[1],0,0,v[0],0,0,v[1],0,
		0,0,-v[0]-v[1],0,0,v[0],0,0,v[1];

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> MuAiT = mStiffness*mArea*Ai.transpose();

	int idx[3] = {mi0,mi1,mi2};

	for(int i=0;i<3;i++)
	{
		J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+0, 3*(index)+0, MuAiT(3*i+0,0)));
		J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+0, 3*(index)+1, MuAiT(3*i+0,1)));
		J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+0, 3*(index)+2, MuAiT(3*i+0,2)));
		J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+1, 3*(index)+0, MuAiT(3*i+1,0)));
		J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+1, 3*(index)+1, MuAiT(3*i+1,1)));
		J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+1, 3*(index)+2, MuAiT(3*i+1,2)));
		J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+2, 3*(index)+0, MuAiT(3*i+2,0)));
		J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+2, 3*(index)+1, MuAiT(3*i+2,1)));
		J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+2, 3*(index)+2, MuAiT(3*i+2,2)));
	}
}

template <typename TinyScalar, typename TinyConstants> 
void TriangleMuscleConstraint<TinyScalar, TinyConstants>::
EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets)
{
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> Ai(3,9);

	Eigen::Matrix<TinyScalar, 2, 1> v = mInvDm*mFiberDirection;
	Ai<<
		-v[0]-v[1],0,0,v[0],0,0,v[1],0,0,
		0,-v[0]-v[1],0,0,v[0],0,0,v[1],0,
		0,0,-v[0]-v[1],0,0,v[0],0,0,v[1];

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> MuAiTAi = mStiffness*mArea*((Ai.transpose())*Ai);

	int idx[3] = {mi0,mi1,mi2};
	for(int i =0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
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
void TriangleMuscleConstraint<TinyScalar, TinyConstants>::
EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x)
{
	Eigen::Matrix<TinyScalar, 3, 1> x0(x.segment<3>(mi0*3));

	Eigen::Matrix32d Ds, P;
	Ds.col(0) = x.segment<3>(mi1*3) - x0;
	Ds.col(1) = x.segment<3>(mi2*3) - x0;

	P.col(0) = Ds.col(0).normalized();
	P.col(1) = (Ds.col(1)-Ds.col(1).dot(P.col(0))*P.col(0)).normalized();

	Eigen::Matrix2d F = P.transpose()*Ds*mInvDm;

	md = (1.0-mActivationLevel)*P*F*mFiberDirection;	
}

template <typename TinyScalar, typename TinyConstants> 
void TriangleMuscleConstraint<TinyScalar, TinyConstants>::
GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d)
{
	d.segment<3>(3*(index)) = md;
	index++;
}

};


#endif