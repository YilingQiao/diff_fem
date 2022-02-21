#ifndef __NORMALIZER_H__
#define __NORMALIZER_H__
#include <Eigen/Core>
#include <vector>

template <typename TinyScalar, typename TinyConstants> 
class Normalizer
{
public: 
	Normalizer(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& real_val_max,const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& real_val_min,
		const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& norm_val_max,const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& norm_val_min);

public:
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> RealToNorm(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& val);
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> NormToReal(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& val);

private:
	int mDim;

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> mRealValMax;	
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> mRealValMin;
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> mRealValDiff;
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> mRealValDiffInv;

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> mNormValMax;	
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> mNormValMin;
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> mNormValDiff;
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> mNormValDiffInv;
};






#include <iostream>
template <typename TinyScalar, typename TinyConstants> 
Normalizer<TinyScalar, TinyConstants>::
Normalizer(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& real_val_max,const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& real_val_min,
		const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& norm_val_max,const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& norm_val_min)
{
	mDim = real_val_max.size();

	mRealValMax.resize(mDim);
	mRealValMin.resize(mDim);

	mNormValMax.resize(mDim);
	mNormValMin.resize(mDim);

	mRealValDiff.resize(mDim);
	mNormValDiff.resize(mDim);

	mRealValDiffInv.resize(mDim);
	mNormValDiffInv.resize(mDim);

	mRealValMax = real_val_max;
	mRealValMin	= real_val_min;

	mNormValMax = norm_val_max;
	mNormValMin	= norm_val_min;

	mRealValDiff = mRealValMax-mRealValMin;
	mNormValDiff = mNormValMax-mNormValMin;

	for(int i=0; i<mDim; i++) {
		mRealValDiffInv[i] = 1.0/mRealValDiff[i];
		mNormValDiffInv[i] = 1.0/mNormValDiff[i];
	}
}

template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> Normalizer<TinyScalar, TinyConstants>::
RealToNorm(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& val)
{
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> val_0_1(mDim);
	
	for(int i=0; i<mDim; i++) 
	{	
		val_0_1[i] = (val[i] - mRealValMin[i]) * mRealValDiffInv[i];
		val_0_1[i] = std::min(std::max(val_0_1[i],0.0),1.0);
	}

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> result(mDim);
	for(int i=0; i<mDim; i++) 
	{
		result[i] = mNormValMin[i] + mNormValDiff[i]*val_0_1[i];
	}

	return result;
}

template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> Normalizer<TinyScalar, TinyConstants>::
NormToReal(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& val)
{
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> val_0_1(mDim);
	
	for(int i=0; i<mDim; i++) 
	{	
		val_0_1[i] = (val[i] - mNormValMin[i]) * mNormValDiffInv[i];
		val_0_1[i] = std::min(
			std::max(val_0_1[i],TinyConstants::zero()),
			TinyConstants::one());
	}

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> result(mDim);
	for(int i=0; i<mDim; i++) 
	{
		result[i] = mRealValMin[i] + mRealValDiff[i]*val_0_1[i];
	}

	return result;
}

#endif