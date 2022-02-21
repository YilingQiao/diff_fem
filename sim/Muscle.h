#ifndef __MUSCLE_H__
#define __MUSCLE_H__
#include "MuscleSegment.h"
#include "fem/Mesh/MeshHeaders.h"

template <typename TinyScalar, typename TinyConstants> 
class Muscle
{
public:
	Muscle(int num_sampling,std::vector<Eigen::Matrix<TinyScalar, 3, 1>> point_list);

	void Initialize(FEM::Mesh<TinyScalar, TinyConstants>* mesh,TinyScalar muscle_stiffness);
	void Reset();

	const std::vector<MuscleSegment<TinyScalar, TinyConstants>*>& GetSegments() {return mSegments;};

	const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& GetStarts() {return mStarts;};	
	const std::vector<Eigen::Matrix<TinyScalar, 3, 1>>& GetEnds() {return mEnds;};	

	const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& GetActions() {return mActions;};
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> GetActionUpperBound() {return mActionUpperBound;};
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> GetActionLowerBound() {return mActionLowerBound;};

	void SetKey(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& key);
	void SetActivationLevels(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& action,const int& phase);
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> GetActivationLevels(){return mActivationLevels;};

// private:	
	int 							mNumSampling;
	std::vector<MuscleSegment<TinyScalar, TinyConstants>*> 	mSegments;

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> 				mActivationLevels;

	std::vector<Eigen::Matrix<TinyScalar, 3, 1>> 	mStarts;
	std::vector<Eigen::Matrix<TinyScalar, 3, 1>> 	mEnds;

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>					mActions;
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> 				mActionUpperBound;	
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> 				mActionLowerBound;

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>					mKey;


	bool isInTetra(const Eigen::Matrix<TinyScalar, 3, 1>& p0,const Eigen::Matrix<TinyScalar, 3, 1>& p1,const Eigen::Matrix<TinyScalar, 3, 1>& p2,
		const Eigen::Matrix<TinyScalar, 3, 1>& p3,const Eigen::Matrix<TinyScalar, 3, 1>& start,const Eigen::Matrix<TinyScalar, 3, 1>& end);
	TinyScalar GetSignedVolume(const Eigen::Matrix<TinyScalar, 3, 1>& p0,const Eigen::Matrix<TinyScalar, 3, 1>& p1,
		const Eigen::Matrix<TinyScalar, 3, 1>& p2,const Eigen::Matrix<TinyScalar, 3, 1>& p3);
	bool isLineTriangleIntersect(const Eigen::Matrix<TinyScalar, 3, 1>& p0,const Eigen::Matrix<TinyScalar, 3, 1>& p1,
		const Eigen::Matrix<TinyScalar, 3, 1>& p2,const Eigen::Matrix<TinyScalar, 3, 1>& start,const Eigen::Matrix<TinyScalar, 3, 1>& end);
	Eigen::Matrix<TinyScalar, 3, 3> GetDm(const Eigen::Matrix<TinyScalar, 3, 1>& p0,const Eigen::Matrix<TinyScalar, 3, 1>& p1,const Eigen::Matrix<TinyScalar, 3, 1>& p2,
		const Eigen::Matrix<TinyScalar, 3, 1>& p3);


};




template <typename TinyScalar, typename TinyConstants> 
Muscle<TinyScalar, TinyConstants>::
Muscle(int num_sampling,std::vector<Eigen::Matrix<TinyScalar, 3, 1>> point_list)
	:mNumSampling(num_sampling)
{
	mActions.resize(4);
	mActions.setZero();

	mActionUpperBound.resize(4);
	mActionLowerBound.resize(4);

	mActionUpperBound[0] = 0.25;
	mActionLowerBound[0] = -0.25;

	mActionUpperBound[1] = 1.0;
	mActionLowerBound[1] = 0.0;
	
	mActionUpperBound[2] = 1.0;
	mActionLowerBound[2] = 0.0;

	mActionUpperBound[3] = 10.0;
	mActionLowerBound[3] = 0.0;

	TinyScalar arc_length =0.0;
	for(int i=0; i<point_list.size()-1; i++) {
		arc_length += (point_list[i+1]-point_list[i]).norm();
	}		

	TinyScalar sampling_length = arc_length/mNumSampling;
	TinyScalar tmp_length = 0.0;

	int start_idx = 0;		

	for(int i=0; i<point_list.size(); i++) {
		tmp_length += (point_list[i+1]-point_list[i]).norm();
		if(tmp_length > sampling_length) {
			mStarts.push_back(point_list[start_idx]);
			mEnds.push_back(point_list[i]);
			tmp_length = 0.0;
			start_idx = i;
		}
	}

	mEnds.pop_back();
	mEnds.push_back(point_list[point_list.size()-1]);

	mActivationLevels.resize(mNumSampling);
	mActivationLevels.setZero();
}

template <typename TinyScalar, typename TinyConstants> 
void Muscle<TinyScalar, TinyConstants>::
Initialize(FEM::Mesh<TinyScalar, TinyConstants>* mesh,TinyScalar muscle_stiffness)
{
	const auto& vertices = mesh->GetVertices();
	const auto& tetrahedron = mesh->GetTetrahedrons();

	for(int i=0; i<mNumSampling; i++) 
	{
		mSegments.push_back(new MuscleSegment<TinyScalar, TinyConstants>());
		auto segment = mSegments.back();

		Eigen::Matrix<TinyScalar, 3, 1> start = mStarts[i]; 
		Eigen::Matrix<TinyScalar, 3, 1> end = mEnds[i]; 

		TinyScalar start_min = 1E5;
		TinyScalar end_min = 1E5;
		
		int start_idx = -1,end_idx = -1;

		for(int j=0;j<tetrahedron.size();j++)
		{
			int i0 = tetrahedron[j][0];
			int i1 = tetrahedron[j][1];
			int i2 = tetrahedron[j][2];
			int i3 = tetrahedron[j][3];

			Eigen::Matrix<TinyScalar, 3, 1> p0 = vertices[i0],p1 = vertices[i1],p2 = vertices[i2],p3 = vertices[i3];

			Eigen::Matrix<TinyScalar, 3, 3> Dm;	
			Dm.template block<3,1>(0,0) = p1 - p0;
			Dm.template block<3,1>(0,1) = p2 - p0;
			Dm.template block<3,1>(0,2) = p3 - p0;

			bool is_under_zero = false;
			if(Dm.determinant()<0.0)
			{
				is_under_zero = true;
				i1 = tetrahedron[j][2];
				i2 = tetrahedron[j][1];
				p1 = vertices[i1];
				p2 = vertices[i2];
				Dm.template block<3,1>(0,0) = p1 - p0;
				Dm.template block<3,1>(0,1) = p2 - p0;
				Dm.template block<3,1>(0,2) = p3 - p0;
			}

			Eigen::Matrix<TinyScalar, 3, 1> fiber_direction = (end-start).normalized();
			Eigen::Matrix<TinyScalar, 3, 1> center = (p0+p1+p2+p3)/4.0;
			TinyScalar dist;
			Eigen::Matrix<TinyScalar, 3, 1> v = center-start;

			TinyScalar t = v.dot(fiber_direction);
			Eigen::Matrix<TinyScalar, 3, 1> p = start + t*fiber_direction; // projection of tet center on fibre
			if(std::max((p-start).norm(),(p-end).norm())/(end-start).norm() < 1.0) // on this fibre
				dist = (p-center).norm();
			else
				dist = 1E5;

			dist = std::min(dist,std::min((center-start).norm(),(center-end).norm())); // dist from center to fibre
			TinyScalar weight = TinyConstants::exp1(-TinyConstants::abs1(dist)/0.005);
	
			if(isInTetra(p0,p1,p2,p3,start,end)) {
				segment->AddMuscleConstraint(new FEM::LinearMuscleConstraint<TinyScalar, TinyConstants>(muscle_stiffness,fiber_direction,0.0,i0,i1,i2,i3,1.0/6.0*Dm.determinant(),Dm.inverse(),weight));
			}
			else if(weight>0.2) {
				segment->AddMuscleConstraint(new FEM::LinearMuscleConstraint<TinyScalar, TinyConstants>(muscle_stiffness,fiber_direction,0.0,i0,i1,i2,i3,1.0/6.0*Dm.determinant(),Dm.inverse(),weight));
			}

			TinyScalar D = Dm.determinant();

			TinyScalar S0,S1,S2,S3;
			S0 = GetDm(start,p1,p2,p3).determinant()/D;
			S1 = GetDm(start,p2,p3,p0).determinant()/D;
			S2 = GetDm(start,p3,p0,p1).determinant()/D;
			S3 = GetDm(start,p0,p1,p2).determinant()/D;

			TinyScalar E0,E1,E2,E3;
			E0 = GetDm(end,p1,p2,p3).determinant()/D;
			E1 = GetDm(end,p2,p3,p0).determinant()/D;
			E2 = GetDm(end,p3,p0,p1).determinant()/D;
			E3 = GetDm(end,p0,p1,p2).determinant()/D;

			if( TinyConstants::abs1(S0+S1+S2+S3-1.0)<1E-4)
			{
				start_idx = j;
				segment->SetStart(Eigen::Vector4i(i0,i1,i2,i3),Eigen::Matrix<TinyScalar, 4, 1>(S0,S1,S2,S3));
				segment->AddMuscleConstraint(new FEM::LinearMuscleConstraint<TinyScalar, TinyConstants>(muscle_stiffness,fiber_direction,0.0,i0,i1,i2,i3,1.0/6.0*Dm.determinant(),Dm.inverse(),weight));
			}

			if( TinyConstants::abs1(E0+E1+E2+E3-1.0)<1E-4)
			{
				end_idx = j;
				segment->SetEnd(Eigen::Vector4i(i0,i1,i2,i3),Eigen::Matrix<TinyScalar, 4, 1>(E0,E1,E2,E3));
				segment->AddMuscleConstraint(new FEM::LinearMuscleConstraint<TinyScalar, TinyConstants>(muscle_stiffness,fiber_direction,0.0,i0,i1,i2,i3,1.0/6.0*Dm.determinant(),Dm.inverse(),weight));
			}
		}
		if(segment->GetStartBarycentric().norm()<1E-6)
		{
			segment->SetStart(segment->GetEndIdx(),segment->GetEndBarycentric());
		}
		else if(segment->GetEndBarycentric().norm()<1E-6)
		{
			segment->SetEnd(segment->GetStartIdx(),segment->GetStartBarycentric());
		}
	}
}

template <typename TinyScalar, typename TinyConstants> 
void Muscle<TinyScalar, TinyConstants>::
SetKey(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& key)
{
	mKey=key;
}

template <typename TinyScalar, typename TinyConstants> 
void Muscle<TinyScalar, TinyConstants>::
SetActivationLevels(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& action,const int& phase)
{
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> new_act(mSegments.size());
	new_act.setZero();

	TinyScalar key = mKey[phase%mKey.size()];

	TinyScalar delta_act = action[0];
	TinyScalar alpha = action[1];
	TinyScalar beta = 1.0-action[1];
	TinyScalar prop_speed = action[3];

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> prev_act = mActivationLevels;

	for(int i=0; i<prop_speed; i++) 
	{
		for(int j=0; j<mNumSampling; j++) {
			if(j == 0)  
				new_act[j] = key + delta_act;
			else
				new_act[j] = alpha*prev_act[j-1] + beta*prev_act[j];
		} 

		prev_act = new_act;
	}

	for(int i=0; i<new_act.size(); i++)
		new_act[i] = std::min(
			std::max(new_act[i],TinyConstants::zero()),
			TinyConstants::one());

	mActivationLevels = new_act;

	for(int i=0;i<mSegments.size();i++)
		mSegments[i]->SetActivationLevel(mActivationLevels[i]);
}

template <typename TinyScalar, typename TinyConstants> 
void Muscle<TinyScalar, TinyConstants>::
Reset()
{
	mActivationLevels.setZero();
	
	for(int i =0;i<mNumSampling;i++)
	{
		mSegments[i]->SetActivationLevel(mActivationLevels[i]);
	}
}

template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, 3, 3> Muscle<TinyScalar, TinyConstants>::
GetDm(const Eigen::Matrix<TinyScalar, 3, 1>& p0,const Eigen::Matrix<TinyScalar, 3, 1>& p1,const Eigen::Matrix<TinyScalar, 3, 1>& p2,const Eigen::Matrix<TinyScalar, 3, 1>& p3)
{
	Eigen::Matrix<TinyScalar, 3, 3> Dm;

	Dm.template block<3,1>(0,0) = p1 - p0;
	Dm.template block<3,1>(0,1) = p2 - p0;
	Dm.template block<3,1>(0,2) = p3 - p0;

	if(Dm.determinant()<0.0)
	{
		Dm.template block<3,1>(0,0) = p2 - p0;
		Dm.template block<3,1>(0,1) = p1 - p0;
		Dm.template block<3,1>(0,2) = p3 - p0;
	}	

	return Dm;
}

template <typename TinyScalar, typename TinyConstants> 
bool Muscle<TinyScalar, TinyConstants>::
isInTetra(const Eigen::Matrix<TinyScalar, 3, 1>& p0,const Eigen::Matrix<TinyScalar, 3, 1>& p1,const Eigen::Matrix<TinyScalar, 3, 1>& p2,const Eigen::Matrix<TinyScalar, 3, 1>& p3,const Eigen::Matrix<TinyScalar, 3, 1>& start,const Eigen::Matrix<TinyScalar, 3, 1>& end)
{
	if(isLineTriangleIntersect(p0,p1,p2,start,end)) {
		return true;
	}
	if(isLineTriangleIntersect(p0,p1,p3,start,end)) {
		return true;
	}
	if(isLineTriangleIntersect(p0,p2,p3,start,end)) {
		return true;	
	}
	if(isLineTriangleIntersect(p1,p2,p3,start,end)) {
		return true;
	}

	return false;
}

template <typename TinyScalar, typename TinyConstants> 
TinyScalar Muscle<TinyScalar, TinyConstants>::
GetSignedVolume(const Eigen::Matrix<TinyScalar, 3, 1>& p0,const Eigen::Matrix<TinyScalar, 3, 1>& p1,const Eigen::Matrix<TinyScalar, 3, 1>& p2,const Eigen::Matrix<TinyScalar, 3, 1>& p3)
{
	return ((p1-p0).cross(p2-p0)).dot(p3-p0)/6.0;
}

template <typename TinyScalar, typename TinyConstants> 
bool Muscle<TinyScalar, TinyConstants>::
isLineTriangleIntersect(const Eigen::Matrix<TinyScalar, 3, 1>& p0,const Eigen::Matrix<TinyScalar, 3, 1>& p1,const Eigen::Matrix<TinyScalar, 3, 1>& p2,const Eigen::Matrix<TinyScalar, 3, 1>& start,const Eigen::Matrix<TinyScalar, 3, 1>& end)
{
	if(GetSignedVolume(start,p0,p1,p2)*GetSignedVolume(end,p0,p1,p2) < 0.0) {
		if(GetSignedVolume(start,end,p0,p1)>0.0 && GetSignedVolume(start,end,p1,p2)>0.0 && GetSignedVolume(start,end,p2,p0)>0.0) {
			return true;
		}
		if(GetSignedVolume(start,end,p0,p1)<0.0 && GetSignedVolume(start,end,p1,p2)<0.0 && GetSignedVolume(start,end,p2,p0)<0.0) {
			return true;
		}
	}

	return false;
}

#endif