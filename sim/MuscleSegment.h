#ifndef __MUSCLE_SEGMENT_H__
#define __MUSCLE_SEGMENT_H__
#include "fem/Constraint/ConstraintHeader.h"
template <typename TinyScalar, typename TinyConstants> 
class MuscleSegment
{
public:
	MuscleSegment();
	void AddMuscleConstraint(FEM::LinearMuscleConstraint<TinyScalar, TinyConstants>* lmc);
	const std::vector<FEM::LinearMuscleConstraint<TinyScalar, TinyConstants>*>& GetMuscleConstraints(){return mMuscleConstraints;};

	void SetActivationLevel(TinyScalar act);
	TinyScalar GetActivationLevel(){return mActivationLevel;};

	void SetStart(const Eigen::Vector4i& start,const Eigen::Matrix<TinyScalar, 4, 1>& barycentric);
	void SetEnd(const Eigen::Vector4i& end,const Eigen::Matrix<TinyScalar, 4, 1>& barycentric);
	
	const Eigen::Vector4i& GetStartIdx() {return mStart;};
	const Eigen::Vector4i& GetEndIdx() {return mEnd;};
	const Eigen::Matrix<TinyScalar, 4, 1>& GetStartBarycentric() {return mStartBarycentric;};
	const Eigen::Matrix<TinyScalar, 4, 1>& GetEndBarycentric() {return mEndBarycentric;};

private:
	TinyScalar 										mActivationLevel;
	std::vector<FEM::LinearMuscleConstraint<TinyScalar, TinyConstants>*> 	mMuscleConstraints;

	Eigen::Vector4i mStart;
	Eigen::Vector4i mEnd;
	Eigen::Matrix<TinyScalar, 4, 1> mStartBarycentric;
	Eigen::Matrix<TinyScalar, 4, 1> mEndBarycentric;
};




template <typename TinyScalar, typename TinyConstants> 
MuscleSegment<TinyScalar, TinyConstants>::
MuscleSegment()
	:mActivationLevel(0.0)
{
	mStart.setZero();
	mStartBarycentric.setZero();
	mEnd.setZero();
	mEndBarycentric.setZero();
}

template <typename TinyScalar, typename TinyConstants> 
void MuscleSegment<TinyScalar, TinyConstants>::
AddMuscleConstraint(FEM::LinearMuscleConstraint<TinyScalar, TinyConstants>* lmc)
{
	mMuscleConstraints.push_back(lmc);	
}

template <typename TinyScalar, typename TinyConstants> 
void MuscleSegment<TinyScalar, TinyConstants>::
SetStart(const Eigen::Vector4i& start,const Eigen::Matrix<TinyScalar, 4, 1>& barycentric)
{
	mStart = start;
	mStartBarycentric = barycentric;
}

template <typename TinyScalar, typename TinyConstants> 
void MuscleSegment<TinyScalar, TinyConstants>::
SetEnd(const Eigen::Vector4i& end,const Eigen::Matrix<TinyScalar, 4, 1>& barycentric)
{
	mEnd = end;
	mEndBarycentric = barycentric;
}

template <typename TinyScalar, typename TinyConstants> 
void MuscleSegment<TinyScalar, TinyConstants>::
SetActivationLevel(TinyScalar activation)
{
	mActivationLevel = activation;
	for(const auto& c : mMuscleConstraints)
	{
		c->SetActivationLevel(activation);
	}
}


#endif