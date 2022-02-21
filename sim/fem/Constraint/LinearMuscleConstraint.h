#ifndef __LINEAR_MUSCLE_CONSTRAINT_H__
#define	__LINEAR_MUSCLE_CONSTRAINT_H__

#include "Constraint.h"
#include "Tensor.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <iostream>

namespace FEM
{
template <typename TinyScalar, typename TinyConstants> 
class LinearMuscleConstraint : public Constraint<TinyScalar, TinyConstants>
{
public:
	LinearMuscleConstraint(const TinyScalar& stiffness,
		const Eigen::Matrix<TinyScalar, 3, 1>& fiber_direction,
		const TinyScalar& activation_level, 
		int i0,int i1,int i2,int i3,
		TinyScalar vol,const Eigen::Matrix<TinyScalar, 3, 3>& invDm,
		// const Eigen::Matrix<TinyScalar, 4, 1>& barycentric1,const Eigen::Matrix<TinyScalar, 4, 1>& barycentric2,
		TinyScalar weight);

protected:
	TinyScalar				mStiffness;
	Eigen::Matrix<TinyScalar, 3, 1>		mFiberDirection;
	TinyScalar 				mActivationLevel;
	int 				mActivationIndex;
	Eigen::Matrix<TinyScalar, 3, 3> 	mddT;

	int 				mi0,mi1,mi2,mi3;
	TinyScalar 				mVol;
	Eigen::Matrix<TinyScalar, 3, 3> 	mInvDm;
	Eigen::Matrix<TinyScalar, 3, 3> 	mDs;
	Eigen::Matrix<TinyScalar, 3, 3> 	mF;

	Eigen::Matrix<TinyScalar, 3, 1>		mp0;

	Eigen::Matrix<TinyScalar, 3, 3>		md;

	TinyScalar				mWeight;
	std::vector<Eigen::Matrix<TinyScalar, 4, 1>> mBarycentric;

public:
	int GetI0() {return mi0;}
	int GetI1() {return mi1;}
	int GetI2() {return mi2;}
	int GetI3() {return mi3;}
	
	int GetDof() override;
	ConstraintType GetType() override;

	void	EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets) override;
	void	EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets) override;
	void 	EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
	void 	GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d);

	void	fixIndex(int offset);
private:
	void	ComputeF(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
	void	ComputeP(Eigen::Matrix<TinyScalar, 3, 3>& P);
	void	ComputedPdF(Tensor3333& dPdF);

	void	Computep0();

public:
	void SetActivationLevel(const TinyScalar& a);
	const TinyScalar& GetActivationLevel();
	const Eigen::Matrix<TinyScalar, 3, 1>& GetFiberDirection();
	int GetActivationIndex(){return mActivationIndex;};
	void SetActivationIndex(const int& i);
	std::vector<Eigen::Matrix<TinyScalar, 4, 1>> GetBarycentric() {return mBarycentric;};
	TinyScalar GetWeight() {return mWeight;};
};




template <typename TinyScalar, typename TinyConstants> 
void LinearMuscleConstraint<TinyScalar, TinyConstants>::fixIndex(int offset) {
  mi0 += offset;
  mi1 += offset;
  mi2 += offset;
  mi3 += offset;
}
template <typename TinyScalar, typename TinyConstants> 
LinearMuscleConstraint<TinyScalar, TinyConstants>::
LinearMuscleConstraint(const TinyScalar& stiffness,
		const Eigen::Matrix<TinyScalar, 3, 1>& fiber_direction,
		const TinyScalar& activation_level, 
		int i0,int i1,int i2,int i3,
		TinyScalar vol,const Eigen::Matrix<TinyScalar, 3, 3>& invDm,
		// const Eigen::Matrix<TinyScalar, 4, 1>& barycentric1,const Eigen::Matrix<TinyScalar, 4, 1>& barycentric2,
		TinyScalar weight)
	:Constraint<TinyScalar, TinyConstants>(stiffness),
	mi0(i0),mi1(i1),mi2(i2),mi3(i3),
	mStiffness(stiffness),mFiberDirection(fiber_direction),
	mVol(vol),mInvDm(invDm),mDs(Eigen::Matrix<TinyScalar, 3, 3>::Zero()),
	mActivationLevel(activation_level),mWeight(weight)
{
	mF.setZero();
	// mBarycentric.clear();
	// mBarycentric.push_back(barycentric1);
	// mBarycentric.push_back(barycentric2);

	mStiffness *= mWeight;
}

template <typename TinyScalar, typename TinyConstants> 
void LinearMuscleConstraint<TinyScalar, TinyConstants>::
ComputeF
(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x)
{
	Eigen::Matrix<TinyScalar, 3, 1> x0(x.template block<3,1>(mi0*3,0));

	Eigen::Matrix<TinyScalar, 3, 3> Ds;

	Ds.template block<3,1>(0,0) = x.template block<3,1>(mi1*3,0)-x0;
	Ds.template block<3,1>(0,1) = x.template block<3,1>(mi2*3,0)-x0;
	Ds.template block<3,1>(0,2) = x.template block<3,1>(mi3*3,0)-x0;

	mDs = Ds;
	mF = mDs * mInvDm;
}

template <typename TinyScalar, typename TinyConstants> 
void LinearMuscleConstraint<TinyScalar, TinyConstants>::
ComputeP
(Eigen::Matrix<TinyScalar, 3, 3>& P)
{
	P = mStiffness*(mF*mFiberDirection*mFiberDirection.transpose() - mp0*mFiberDirection.transpose());
}	

template <typename TinyScalar, typename TinyConstants> 
void LinearMuscleConstraint<TinyScalar, TinyConstants>::
ComputedPdF
(Tensor3333& dPdF)
{
	Tensor3333 dFdF;
	dFdF.SetIdentity();

	for(int i=0; i<3;i++)
		for(int j=0;j<3;j++)
			dPdF(i,j) = mStiffness*dFdF(i,j)*mFiberDirection*mFiberDirection.transpose();
}

template <typename TinyScalar, typename TinyConstants> 
void LinearMuscleConstraint<TinyScalar, TinyConstants>::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets)
{
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> Ai(3,12);

	Eigen::Matrix<TinyScalar, 3, 1> v = mInvDm*mFiberDirection;

	TinyScalar a,b,c;
	a = v[0];
	b = v[1];
	c = v[2];

	Ai<<
		-(a+b+c),0,0,a,0,0,b,0,0,c,0,0,
		0,-(a+b+c),0,0,a,0,0,b,0,0,c,0,
		0,0,-(a+b+c),0,0,a,0,0,b,0,0,c;

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> MuAiT = mVol*mStiffness*Ai.transpose();

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
void LinearMuscleConstraint<TinyScalar, TinyConstants>::	
EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets)
{
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> Ai(3,12);
	Eigen::Matrix<TinyScalar, 3, 1> v = mInvDm*mFiberDirection;

	TinyScalar a,b,c;
	a = v[0];
	b = v[1];
	c = v[2];

	Ai<<
		-(a+b+c),0,0,a,0,0,b,0,0,c,0,0,
		0,-(a+b+c),0,0,a,0,0,b,0,0,c,0,
		0,0,-(a+b+c),0,0,a,0,0,b,0,0,c;

	auto MuAiTAi = mVol*mStiffness*((Ai.transpose())*Ai);

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
void LinearMuscleConstraint<TinyScalar, TinyConstants>::	
EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x)
{
	ComputeF(x);
	Computep0();
}

template <typename TinyScalar, typename TinyConstants> 
void LinearMuscleConstraint<TinyScalar, TinyConstants>::
GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d)
{
	d.template block<3,1>(3*index,0) = mp0;
	index++;
}

template <typename TinyScalar, typename TinyConstants> 
void LinearMuscleConstraint<TinyScalar, TinyConstants>::
Computep0()
{
	mp0 = (1.0-mActivationLevel)*mF*mFiberDirection;
	// mp0 = (1.0)*mF*mFiberDirection;
}	

template <typename TinyScalar, typename TinyConstants> 
int LinearMuscleConstraint<TinyScalar, TinyConstants>::
GetDof()
{
	return 1;
}

template <typename TinyScalar, typename TinyConstants> 
ConstraintType LinearMuscleConstraint<TinyScalar, TinyConstants>::
GetType()
{
	return ConstraintType::LINEAR_MUSCLE;
}

template <typename TinyScalar, typename TinyConstants> 
void LinearMuscleConstraint<TinyScalar, TinyConstants>::
SetActivationLevel(const TinyScalar& a) 
{
	mActivationLevel = a;
} 

template <typename TinyScalar, typename TinyConstants> 
const TinyScalar& LinearMuscleConstraint<TinyScalar, TinyConstants>::
GetActivationLevel() 
{
	return mActivationLevel;
}

template <typename TinyScalar, typename TinyConstants> 
const Eigen::Matrix<TinyScalar, 3, 1>& LinearMuscleConstraint<TinyScalar, TinyConstants>::
GetFiberDirection()
{
	return mFiberDirection;
}

template <typename TinyScalar, typename TinyConstants> 
void LinearMuscleConstraint<TinyScalar, TinyConstants>::
SetActivationIndex(const int& i)
{
	mActivationIndex = i;
}

};

#endif