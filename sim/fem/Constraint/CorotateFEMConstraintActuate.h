#ifndef __COROTATE_FEM_CONSTRAINT_ACTUATE_H__
#define	__COROTATE_FEM_CONSTRAINT_ACTUATE_H__
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
class CorotateFEMConstraintActuate : public Constraint<TinyScalar, TinyConstants>
{
public:
	CorotateFEMConstraintActuate(
		const TinyScalar& stiffness,
		const TinyScalar& poisson_ratio,
		int i0,int i1,int i2,int i3,
		TinyScalar volume,const Eigen::Matrix<TinyScalar, 3, 3>& invDm);

	int GetI0() {return mi0;}
	int GetI1() {return mi1;}
	int GetI2() {return mi2;}
	int GetI3() {return mi3;}

	int GetDof() override;
	ConstraintType GetType() override;

	void SetActivationLevel(const TinyScalar& a);

private:
	void	ComputeF(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
	void	ComputeP(Eigen::Matrix<TinyScalar, 3, 3>& P);
	void	ComputedPdF(Tensor3333& dPdF);
	void	ComputeSVD(const Eigen::Matrix<TinyScalar, 3, 3>& F);

	void	EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets);
	void	EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets);
	void 	EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
	void 	GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d);

protected:
	int 				mi0,mi1,mi2,mi3;
	TinyScalar 				mVol;
	TinyScalar 				mMu,mLambda;
	TinyScalar 				mPoissonRatio;
	Eigen::Matrix<TinyScalar, 3, 3> 	mInvDm;
	Eigen::Matrix<TinyScalar, 3, 3> 	mDs; // [x1-x0; x2-x0; x3-x0]
	Eigen::Matrix<TinyScalar, 3, 3> 	mF;
	Eigen::Matrix<TinyScalar, 3, 3> 	mR,mU,mV,mD; // R=UV;  D=svd.singularValues();

	Eigen::Matrix<TinyScalar, 3, 3>		md;
	Eigen::Matrix<TinyScalar, 3, 3>		md_volume;

	TinyScalar 				mActivationLevel = 1;
};


#define EPS 1E-4

template <typename TinyScalar, typename TinyConstants> 
CorotateFEMConstraintActuate<TinyScalar, TinyConstants>::
CorotateFEMConstraintActuate(const TinyScalar& stiffness,const TinyScalar& poisson_ratio,
	int i0,int i1,int i2,int i3,TinyScalar vol,const Eigen::Matrix<TinyScalar, 3, 3>& invDm)
	:Constraint<TinyScalar, TinyConstants>(stiffness),
	mPoissonRatio(poisson_ratio),
	mi0(i0),mi1(i1),mi2(i2),mi3(i3),
	mMu(stiffness/((1.0+poisson_ratio))),
	mLambda(stiffness*poisson_ratio/((1.0+poisson_ratio)*(1-2.0*poisson_ratio))),
	mVol(vol),mInvDm(invDm),mDs(Eigen::Matrix<TinyScalar, 3, 3>::Zero())
{
	mF.setZero();
	mR.setZero();
	mU.setZero();	
	mV.setZero();
}

template <typename TinyScalar, typename TinyConstants> 
void CorotateFEMConstraintActuate<TinyScalar, TinyConstants>::
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

	ComputeSVD(mF);
}

template <typename TinyScalar, typename TinyConstants> 
void CorotateFEMConstraintActuate<TinyScalar, TinyConstants>::
ComputeP
(Eigen::Matrix<TinyScalar, 3, 3>& P)
{
	P = mMu*(mF - mR)
		+ mLambda*((mR.transpose()*mF-Eigen::Matrix<TinyScalar, 3, 3>::Identity()).trace())*mR;
}

template <typename TinyScalar, typename TinyConstants> 
void CorotateFEMConstraintActuate<TinyScalar, TinyConstants>::
ComputedPdF
(Tensor3333& dPdF)
{
	Tensor3333 dFdF, dRdF;
	dFdF.SetIdentity();

	for(int i =0;i<3;i++) {
		for(int j=0;j<3;j++) {
			Eigen::Matrix<TinyScalar, 3, 3> M = mU.transpose()*dFdF(i,j)*mV;

			if(fabs(mD(0,0)-mD(1,1)) < EPS && fabs(mD(0,0)-mD(2,2)) < EPS) {
				Eigen::Matrix<TinyScalar, 3, 3> off_diag_M;
				off_diag_M.setZero();
				for(int a=0; a<3; a++) {
					for(int b=0; b<3; b++) {
						if(a==b)
							continue;
						else
							off_diag_M(a,b) = M(a,b) / mD(0,0);
					}
				}

				dRdF(i,j) = mU*off_diag_M*mV.transpose();
			} else {
				Eigen::Matrix<TinyScalar, 2, 1> unknown_side, known_side;
				Eigen::Matrix2d known_matrix;
				Eigen::Matrix<TinyScalar, 3, 3> U_tilde, V_tilde;
				U_tilde.setZero(); 
				V_tilde.setZero();
				Eigen::Matrix2d reg;
				reg.setZero();
				reg(0,0) = reg(1,1) = EPS;
				for (unsigned int row = 0; row < 3; row++) {
					for (unsigned int col = 0; col < row; col++) {
						known_side = Eigen::Matrix<TinyScalar, 2, 1>(M(col, row), M(row, col));
						known_matrix.block<2, 1>(0, 0) = Eigen::Matrix<TinyScalar, 2, 1>(-mD(row,row), mD(col,col));
						known_matrix.block<2, 1>(0, 1) = Eigen::Matrix<TinyScalar, 2, 1>(-mD(col,col), mD(row,row));

						if (fabs(mD(row,row) - mD(col,col) < EPS))
							known_matrix += reg;
						else
							assert(fabs(known_matrix.determinant()) > 1E-6);

						unknown_side = known_matrix.inverse() * known_side;
						U_tilde(row, col) = unknown_side[0];
						U_tilde(col, row) = -U_tilde(row, col);
						V_tilde(row, col) = unknown_side[1];
						V_tilde(col, row) = -V_tilde(row, col);
					}
				}
				Eigen::Matrix<TinyScalar, 3, 3> deltaU = mU*U_tilde;
				Eigen::Matrix<TinyScalar, 3, 3> deltaV = V_tilde*mV.transpose();

				dRdF(i, j) = deltaU*mV.transpose() + mU*deltaV;
			}
		}
	}

	Tensor3333 lambda_term;
	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			lambda_term(i,j) =
				(dRdF(i,j).transpose()*mF+mR.transpose()*dFdF(i,j)).trace()*mR +
				(mR.transpose()*mF-Eigen::Matrix<TinyScalar, 3, 3>::Identity()).trace()*dRdF(i,j);
		}
	}

	dPdF = (dFdF-dRdF)*mMu + mLambda*lambda_term;
}

template <typename TinyScalar, typename TinyConstants> 
void CorotateFEMConstraintActuate<TinyScalar, TinyConstants>::
ComputeSVD(const Eigen::Matrix<TinyScalar, 3, 3>& F)
{
	// #pragma omp critical 
	// {
		Eigen::JacobiSVD<Eigen::Matrix<TinyScalar, 3, 3>> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Matrix<TinyScalar, 3, 1> D = svd.singularValues();
		
		mD.setZero();	

		mD(0,0) = D[0];
		mD(1,1) = D[1];
		mD(2,2) = D[2];

		mU = svd.matrixU();
		mV = svd.matrixV();
		mR = mU*mV.transpose();
		mF = F;
	// }
}

template <typename TinyScalar, typename TinyConstants> 
int CorotateFEMConstraintActuate<TinyScalar, TinyConstants>::
GetDof()
{
	return 6;
}

template <typename TinyScalar, typename TinyConstants> 
ConstraintType CorotateFEMConstraintActuate<TinyScalar, TinyConstants>::
GetType()
{
	return ConstraintType::PNEUMATIC;
}

template <typename TinyScalar, typename TinyConstants> 
void CorotateFEMConstraintActuate<TinyScalar, TinyConstants>::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets)
{
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> Ai(3*3,3*4);
	TinyScalar d11 = mInvDm(0,0);
	TinyScalar d12 = mInvDm(0,1);
	TinyScalar d13 = mInvDm(0,2);
	TinyScalar d21 = mInvDm(1,0);
	TinyScalar d22 = mInvDm(1,1);
	TinyScalar d23 = mInvDm(1,2);
	TinyScalar d31 = mInvDm(2,0);
	TinyScalar d32 = mInvDm(2,1);
	TinyScalar d33 = mInvDm(2,2);

	Ai<<
		-d11-d21-d31,0,0,d11,0,0,d21,0,0,d31,0,0,
		0,-d11-d21-d31,0,0,d11,0,0,d21,0,0,d31,0,
		0,0,-d11-d21-d31,0,0,d11,0,0,d21,0,0,d31,
		-d12-d22-d32,0,0,d12,0,0,d22,0,0,d32,0,0,
		0,-d12-d22-d32,0,0,d12,0,0,d22,0,0,d32,0,
		0,0,-d12-d22-d32,0,0,d12,0,0,d22,0,0,d32,
		-d13-d23-d33,0,0,d13,0,0,d23,0,0,d33,0,0,
		0,-d13-d23-d33,0,0,d13,0,0,d23,0,0,d33,0,
		0,0,-d13-d23-d33,0,0,d13,0,0,d23,0,0,d33;

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> MuAiT = mMu*mVol*Ai.transpose();
	int idx[4] = {mi0,mi1,mi2,mi3};
	for(int i =0;i<4;i++)
	{
		for(int j=0;j<3;j++)
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
	index+=3;

	MuAiT = (MuAiT*mPoissonRatio).eval();
	for(int i =0;i<4;i++)
	{
		for(int j=0;j<3;j++)
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
	index+=3;
}

template <typename TinyScalar, typename TinyConstants> 
void CorotateFEMConstraintActuate<TinyScalar, TinyConstants>::
EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets)
{
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> Ai(3*3,3*4);
	TinyScalar d11 = mInvDm(0,0);
	TinyScalar d12 = mInvDm(0,1);
	TinyScalar d13 = mInvDm(0,2);
	TinyScalar d21 = mInvDm(1,0);
	TinyScalar d22 = mInvDm(1,1);
	TinyScalar d23 = mInvDm(1,2);
	TinyScalar d31 = mInvDm(2,0);
	TinyScalar d32 = mInvDm(2,1);
	TinyScalar d33 = mInvDm(2,2);

	Ai<<
		-d11-d21-d31,0,0,d11,0,0,d21,0,0,d31,0,0,
		0,-d11-d21-d31,0,0,d11,0,0,d21,0,0,d31,0,
		0,0,-d11-d21-d31,0,0,d11,0,0,d21,0,0,d31,
		-d12-d22-d32,0,0,d12,0,0,d22,0,0,d32,0,0,
		0,-d12-d22-d32,0,0,d12,0,0,d22,0,0,d32,0,
		0,0,-d12-d22-d32,0,0,d12,0,0,d22,0,0,d32,
		-d13-d23-d33,0,0,d13,0,0,d23,0,0,d33,0,0,
		0,-d13-d23-d33,0,0,d13,0,0,d23,0,0,d33,0,
		0,0,-d13-d23-d33,0,0,d13,0,0,d23,0,0,d33;

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> MuAiTAi = mMu*mVol*((Ai.transpose())*Ai);
	int idx[4] = {mi0,mi1,mi2,mi3};
	//MuAiT --- 12x12 matrix
	for(int i =0;i<4;i++)
	{
		for(int j=0;j<4;j++)
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

	MuAiTAi = (MuAiTAi*mPoissonRatio).eval();
	for(int i =0;i<4;i++)
	{
		for(int j=0;j<4;j++)
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
void CorotateFEMConstraintActuate<TinyScalar, TinyConstants>::
EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x)
{
	ComputeF(x);
	
	md = mR;
	if(mF.determinant()<0)
		md.template block<3,1>(0,2) = -mR.template block<3,1>(0,2);

	Eigen::Matrix<TinyScalar, 3, 1> S = mD.diagonal();
	Eigen::Matrix<TinyScalar, 3, 1> D;
	D.setZero();
	TinyScalar CD;
	for(int i=0;i<5;i++)
	{
		CD = (S[0]+D[0])*(S[1]+D[1])*(S[2]+D[2])-mActivationLevel;
		Eigen::Matrix<TinyScalar, 3, 1> gradCD( (S[1]+D[1])*(S[2]+D[2]),
								(S[0]+D[0])*(S[2]+D[2]),
								(S[0]+D[0])*(S[1]+D[1]));

		D = (gradCD.dot(D) -CD)/(gradCD.squaredNorm())*gradCD;
	}

	md_volume = mU*((S+D).asDiagonal())*mV.transpose();
}

template <typename TinyScalar, typename TinyConstants> 
void CorotateFEMConstraintActuate<TinyScalar, TinyConstants>::
GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d)
{
	d.template block<3,1>(3*(index+0),0) = md.template block<3,1>(0,0);
	d.template block<3,1>(3*(index+1),0) = md.template block<3,1>(0,1);
	d.template block<3,1>(3*(index+2),0) = md.template block<3,1>(0,2);
	index+=3;

	d.template block<3,1>(3*(index+0),0) = md_volume.template block<3,1>(0,0);
	d.template block<3,1>(3*(index+1),0) = md_volume.template block<3,1>(0,1);
	d.template block<3,1>(3*(index+2),0) = md_volume.template block<3,1>(0,2);
	index+=3;
}


template <typename TinyScalar, typename TinyConstants> 
void CorotateFEMConstraintActuate<TinyScalar, TinyConstants>::
SetActivationLevel(const TinyScalar& a) 
{
	mActivationLevel = a;
}

#undef EPS

};
#endif