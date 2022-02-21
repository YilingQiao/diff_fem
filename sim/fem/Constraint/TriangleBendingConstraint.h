#ifndef __TRIANGLE_BENDING_CONSTRAINT_H__
#define __TRIANGLE_BENDING_CONSTRAINT_H__
#include "Constraint.h"
#include <Eigen/SVD>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <iostream>

namespace Eigen
{
using Matrix34d = Matrix<TinyScalar, 3, 4>;
};
namespace FEM
{
template <typename TinyScalar, typename TinyConstants> 
class TriangleBendingConstraint : public Constraint<TinyScalar, TinyConstants>
{
public:
	TriangleBendingConstraint(TinyScalar stiffness,int i0,int i1,int i2,int i3,TinyScalar voronoi_area,TinyScalar n,const Eigen::Matrix<TinyScalar, 4, 1>& w);
	int GetDof() override;

	void EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets);
	void EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets);
	void EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
	void GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d);

	void	fixIndex(int offset);
	int GetI0() {return mi0;}
	int GetI1() {return mi1;}
	int GetI2() {return mi2;}
	int GetI3() {return mi3;}
protected:
	int mi0,mi1,mi2,mi3;
	TinyScalar mn,mVoronoiArea;
	Eigen::Matrix<TinyScalar, 4, 1> mw;

	Eigen::Matrix<TinyScalar, 3, 1> md;
};



template <typename TinyScalar, typename TinyConstants> 
void TriangleBendingConstraint<TinyScalar, TinyConstants>::fixIndex(int offset) {
  mi0 += offset;
  mi1 += offset;
  mi2 += offset;
  mi3 += offset;
}


template <typename TinyScalar, typename TinyConstants> 
TriangleBendingConstraint<TinyScalar, TinyConstants>::
TriangleBendingConstraint(TinyScalar stiffness,int i0,int i1,int i2,int i3,TinyScalar voronoi_area,TinyScalar n,const Eigen::Matrix<TinyScalar, 4, 1>& w)
	:Constraint<TinyScalar, TinyConstants>(stiffness),mi0(i0),mi1(i1),mi2(i2),mi3(i3),mVoronoiArea(voronoi_area),mn(n),mw(w)
{
}

template <typename TinyScalar, typename TinyConstants> 
int TriangleBendingConstraint<TinyScalar, TinyConstants>::
GetDof()
{
	return 1;
}

template <typename TinyScalar, typename TinyConstants> 
void TriangleBendingConstraint<TinyScalar, TinyConstants>::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets)
{
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> Ai(3,12);

	Ai<<
		mw[0],0,0,mw[1],0,0,mw[2],0,0,mw[3],0,0,
		0,mw[0],0,0,mw[1],0,0,mw[2],0,0,mw[3],0,
		0,0,mw[0],0,0,mw[1],0,0,mw[2],0,0,mw[3];

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> MuAiT = mVoronoiArea*mStiffness*Ai.transpose();

	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+0,3*index+0,MuAiT(3*0+0,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+0,3*index+1,MuAiT(3*0+0,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+0,3*index+2,MuAiT(3*0+0,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+1,3*index+0,MuAiT(3*0+1,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+1,3*index+1,MuAiT(3*0+1,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+1,3*index+2,MuAiT(3*0+1,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+2,3*index+0,MuAiT(3*0+2,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+2,3*index+1,MuAiT(3*0+2,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi0+2,3*index+2,MuAiT(3*0+2,3*0+2)));

	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi1+0,3*index+0,MuAiT(3*1+0,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi1+0,3*index+1,MuAiT(3*1+0,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi1+0,3*index+2,MuAiT(3*1+0,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi1+1,3*index+0,MuAiT(3*1+1,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi1+1,3*index+1,MuAiT(3*1+1,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi1+1,3*index+2,MuAiT(3*1+1,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi1+2,3*index+0,MuAiT(3*1+2,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi1+2,3*index+1,MuAiT(3*1+2,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi1+2,3*index+2,MuAiT(3*1+2,3*0+2)));

	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi2+0,3*index+0,MuAiT(3*2+0,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi2+0,3*index+1,MuAiT(3*2+0,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi2+0,3*index+2,MuAiT(3*2+0,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi2+1,3*index+0,MuAiT(3*2+1,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi2+1,3*index+1,MuAiT(3*2+1,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi2+1,3*index+2,MuAiT(3*2+1,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi2+2,3*index+0,MuAiT(3*2+2,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi2+2,3*index+1,MuAiT(3*2+2,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi2+2,3*index+2,MuAiT(3*2+2,3*0+2)));

	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi3+0,3*index+0,MuAiT(3*3+0,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi3+0,3*index+1,MuAiT(3*3+0,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi3+0,3*index+2,MuAiT(3*3+0,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi3+1,3*index+0,MuAiT(3*3+1,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi3+1,3*index+1,MuAiT(3*3+1,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi3+1,3*index+2,MuAiT(3*3+1,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi3+2,3*index+0,MuAiT(3*3+2,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi3+2,3*index+1,MuAiT(3*3+2,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<TinyScalar>(3*mi3+2,3*index+2,MuAiT(3*3+2,3*0+2)));
}

template <typename TinyScalar, typename TinyConstants> 
void TriangleBendingConstraint<TinyScalar, TinyConstants>::
EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets)
{
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> Ai(3,12);

	Ai<<
		mw[0],0,0,mw[1],0,0,mw[2],0,0,mw[3],0,0,
		0,mw[0],0,0,mw[1],0,0,mw[2],0,0,mw[3],0,
		0,0,mw[0],0,0,mw[1],0,0,mw[2],0,0,mw[3];


	auto MuAiTAi = mVoronoiArea*mStiffness*((Ai.transpose())*Ai);

	int idx[4] = {mi0,mi1,mi2,mi3};

	for(int i =0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+0,3*idx[j]+0,MuAiTAi(3*i+0,3*j+0)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+0,3*idx[j]+1,MuAiTAi(3*i+0,3*j+1)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+0,3*idx[j]+2,MuAiTAi(3*i+0,3*j+2)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+1,3*idx[j]+0,MuAiTAi(3*i+1,3*j+0)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+1,3*idx[j]+1,MuAiTAi(3*i+1,3*j+1)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+1,3*idx[j]+2,MuAiTAi(3*i+1,3*j+2)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+2,3*idx[j]+0,MuAiTAi(3*i+2,3*j+0)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+2,3*idx[j]+1,MuAiTAi(3*i+2,3*j+1)));
			L_triplets.push_back(Eigen::Triplet<TinyScalar>(3*idx[i]+2,3*idx[j]+2,MuAiTAi(3*i+2,3*j+2)));
		}
	}
}

template <typename TinyScalar, typename TinyConstants> 
void TriangleBendingConstraint<TinyScalar, TinyConstants>::
EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x)
{
	Eigen::Matrix<TinyScalar, 3, 1> e = Eigen::Matrix<TinyScalar, 3, 1>::Zero();
	if( mn > 1E-6 )
	{
		e += mw[0] * x.segment<3>(mi0);
		e += mw[1] * x.segment<3>(mi1);
		e += mw[2] * x.segment<3>(mi2);
		e += mw[3] * x.segment<3>(mi3);

		TinyScalar l = e.norm();

		if(l>1E-6)
			e *= mn/l;
	}
	md = e;
}

template <typename TinyScalar, typename TinyConstants> 
void TriangleBendingConstraint<TinyScalar, TinyConstants>::
GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d)
{
	d.segment<3>(3*index) = md;
	index++;
}
};


#endif