#ifndef __TRIANGLE_STRAIN_CONSTRAINT_H__
#define __TRIANGLE_STRAIN_CONSTRAINT_H__
#include "Constraint.h"
#include <Eigen/SVD>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <iostream>

namespace Eigen
{
using Matrix32d = Matrix<TinyScalar, 3, 2>;
};
namespace FEM
{
template <typename TinyScalar, typename TinyConstants> 
class TriangleStrainConstraint : public Constraint<TinyScalar, TinyConstants>
{
public:
	TriangleStrainConstraint(TinyScalar stiffness,int i0,int i1,int i2,TinyScalar area,const Eigen::Matrix2d& invDm);
	int GetDof() override;

	void EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets);
	void EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets);
	void EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
	void GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d);

	int GetI0() {return mi0;}
	int GetI1() {return mi1;}
	int GetI2() {return mi2;}
protected:
	int mi0,mi1,mi2;
	TinyScalar mArea;
	Eigen::Matrix2d mInvDm;

	Eigen::Matrix32d md;
};


template <typename TinyScalar, typename TinyConstants> 
TriangleStrainConstraint<TinyScalar, TinyConstants>::
TriangleStrainConstraint(TinyScalar stiffness,int i0,int i1,int i2,TinyScalar area,const Eigen::Matrix2d& invDm)
	:Constraint(stiffness),mi0(i0),mi1(i1),mi2(i2),mArea(area),mInvDm(invDm)
{

}


template <typename TinyScalar, typename TinyConstants> 
int TriangleStrainConstraint<TinyScalar, TinyConstants>::
GetDof()
{
	return 2;
}


template <typename TinyScalar, typename TinyConstants> 
void TriangleStrainConstraint<TinyScalar, TinyConstants>::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets)
{
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> Ai(3*2,3*3);
	TinyScalar d11 = mInvDm(0,0);
	TinyScalar d12 = mInvDm(0,1);
	TinyScalar d21 = mInvDm(1,0);
	TinyScalar d22 = mInvDm(1,1);

	Ai<<
		-d11-d21,0,0,d11,0,0,d21,0,0,
		0,-d11-d21,0,0,d11,0,0,d21,0,
		0,0,-d11-d21,0,0,d11,0,0,d21,
		-d12-d22,0,0,d12,0,0,d22,0,0,
		0,-d12-d22,0,0,d12,0,0,d22,0,
		0,0,-d12-d22,0,0,d12,0,0,d22;

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> MuAiT = mStiffness*mArea*Ai.transpose();
	int idx[3] = {mi0,mi1,mi2};

	for(int i =0;i<3;i++)
	{
		for(int j=0;j<2;j++)
		{
			//MuAiT.block [i,j] -- 3x3 matrix
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

}

template <typename TinyScalar, typename TinyConstants> 
void TriangleStrainConstraint<TinyScalar, TinyConstants>::
EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets)
{
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> Ai(3*2,3*3);
	TinyScalar d11 = mInvDm(0,0);
	TinyScalar d12 = mInvDm(0,1);
	TinyScalar d21 = mInvDm(1,0);
	TinyScalar d22 = mInvDm(1,1);

	Ai<<
		-d11-d21,0,0,d11,0,0,d21,0,0,
		0,-d11-d21,0,0,d11,0,0,d21,0,
		0,0,-d11-d21,0,0,d11,0,0,d21,
		-d12-d22,0,0,d12,0,0,d22,0,0,
		0,-d12-d22,0,0,d12,0,0,d22,0,
		0,0,-d12-d22,0,0,d12,0,0,d22;

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> MuAiTAi = mStiffness*mArea*((Ai.transpose())*Ai);
	int idx[3] = {mi0,mi1,mi2};

	for(int i =0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			//MuAiT.block [i,j] -- 3x3 matrix
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
void TriangleStrainConstraint<TinyScalar, TinyConstants>::
EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x)
{
	Eigen::Matrix<TinyScalar, 3, 1> x0(x.segment<3>(mi0*3));

	Eigen::Matrix32d Ds, P;
	Ds.col(0) = x.segment<3>(mi1*3) - x0;
	Ds.col(1) = x.segment<3>(mi2*3) - x0;

	P.col(0) = Ds.col(0).normalized();
	P.col(1) = (Ds.col(1)-Ds.col(1).dot(P.col(0))*P.col(0)).normalized();

	Eigen::Matrix2d F = P.transpose()*Ds*mInvDm;

	Eigen::JacobiSVD<Eigen::Matrix2d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);

	md = P*(svd.matrixU())*(svd.matrixV().transpose());
}

template <typename TinyScalar, typename TinyConstants> 
void TriangleStrainConstraint<TinyScalar, TinyConstants>::
GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d)
{
	d.segment<3>(3*(index+0)) = md.col(0);
	d.segment<3>(3*(index+1)) = md.col(1);

	index+=2;
}
};


#endif