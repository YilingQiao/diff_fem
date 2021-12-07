#ifndef __ENVIRONMENT_H__
#define __ENVIRONMENT_H__
#include <deque>
#include "fem/World.h"
#include "Dfobj.h"
#include "Dfcloth.h"
#include "Normalizer.h"
// class Dfobj;
// class Normalizer;

template <typename TinyScalar, typename TinyConstants> 
class Environment
{
public:
	Environment();
    Environment(int total_step, const std::string& config_path, bool do_collision);

	FEM::World<TinyScalar, TinyConstants>* GetSoftWorld(){return mSoftWorld;};
	Dfobj<TinyScalar, TinyConstants>* GetDfobj(){return mDfobj;};

	TinyScalar GetSimulationHz(){return mSimulationHz;};
	TinyScalar GetControlHz(){return mControlHz;};

	void Step();
    void Step(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x,
        const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& v,
        const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& jq,
        const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& jqd,
        const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& tau,
        Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& next_x,
        Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& next_v,
        Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& next_jq,
        Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& next_jqd);
    void SaveObj(std::string filename, const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &x);

	void Reset();

	const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& GetStates();

	void InitializeActions();
	const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& GetActions() {return mActions;};
	void SetActions(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& actions);

	std::map<std::string,TinyScalar> GetRewards();

	bool isEndOfEpisode();

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> GetNormLowerBound() {return mNormLowerBound;};
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> GetNormUpperBound() {return mNormUpperBound;};
	Normalizer<TinyScalar, TinyConstants>*	GetNormalizer() {return mNormalizer;};

	void SetPhase(const int& phase); 
	const int& GetPhase() {return mPhase;};

	Eigen::Matrix<TinyScalar, 3, 1> GetAverageVelocity() {return mAverageVelocity;};
	Eigen::Matrix<TinyScalar, 3, 1> GetTargetVelocity() {return mTargetVelocity;};
	void UpdateRandomTargetVelocity();
    void OutputSurface(std::string filename);

	FEM::World<TinyScalar, TinyConstants>*					mSoftWorld;
    Dfcloth<TinyScalar, TinyConstants>*                       mDfcloth;
    Dfobj<TinyScalar, TinyConstants>*                       mDfobj;

	int 							mSimulationHz;
	int 							mControlHz;

	Normalizer<TinyScalar, TinyConstants>*					mNormalizer;
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> 				mNormLowerBound;
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> 				mNormUpperBound;

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>					mStates;
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>					mActions;

	int 							mPhase, mLastPhase;
    int                             mTotalStep;
	Eigen::Matrix<TinyScalar, 3, 1>					mTargetVelocity;
	Eigen::Matrix<TinyScalar, 3, 1>					mAverageVelocity;
	std::deque<Eigen::Matrix<TinyScalar, 3, 1>>		mAverageVelocityDeque;
};




#include <chrono>
#include <random>
#include <fstream>

template <typename TinyScalar, typename TinyConstants> 
Environment<TinyScalar, TinyConstants>::
Environment()
    :mPhase(0), mTotalStep(-1)
{
	mSimulationHz = 200;
    mControlHz = 30;

    mSoftWorld = new FEM::World<TinyScalar, TinyConstants>(
        1.0/mSimulationHz,
        30,
        0.9999,
        true
    );

    // TODO: initialize log

    // mDfobj = new Dfobj<TinyScalar, TinyConstants>(5E5,2E5,0.4);
    // mDfobj->SetDfobj(std::string(SOFTCON_DIR)+"/data/octopus_rigid.meta");
    // mDfobj->Initialize(mSoftWorld);

    mDfobj = new Dfobj<TinyScalar, TinyConstants>(5E5,2E5,0.4);
    mDfobj->SetDfobj(std::string(SOFTCON_DIR)+"/data/fish/goldfish.json");
    mDfobj->Initialize(mSoftWorld);


    // mDfcloth = new Dfcloth<TinyScalar, TinyConstants>(5E5,2E5,0.4);
    // mDfcloth->SetDfobj(std::string(SOFTCON_DIR)+"/data/json_files/sphere_3_layers.json");
    // mDfcloth->Initialize(mSoftWorld);

    

    mSoftWorld->Initialize();
    mDfobj->SetInitReferenceRotation(mSoftWorld->GetPositions());
    std::cout << "------------ 1\n";
    InitializeActions();
    std::cout << "------------ 2\n";

    mTargetVelocity.setZero();
    std::cout << "------------ 3\n";

    mAverageVelocity.setZero();
    mAverageVelocityDeque.clear();

    UpdateRandomTargetVelocity();
    std::cout << "------------ 4\n";
}


template <typename TinyScalar, typename TinyConstants> 
Environment<TinyScalar, TinyConstants>::
Environment(int total_step, const std::string& config_path, bool do_collision)
    :mPhase(0), mTotalStep(total_step), mLastPhase(-1)
{
    std::cout << mTotalStep << " mTotalStep\n";
    std::cout << config_path << " config_path\n";
    mSimulationHz = 200;
    mControlHz = 30;

    mSoftWorld = new FEM::World<TinyScalar, TinyConstants>(
        1.0/mSimulationHz,
        20,
        0.9999,
        do_collision
    );

    // TODO: initialize log

    // mDfobj = new Dfobj<TinyScalar, TinyConstants>(5E5,2E5,0.4);
    // mDfobj->SetDfobj(std::string(SOFTCON_DIR)+"/data/octopus_rigid.meta");
    // mDfobj->Initialize(mSoftWorld);

    mDfobj = new Dfobj<TinyScalar, TinyConstants>(5E5,2E5,0.4);
    mDfobj->SetDfobj(std::string(SOFTCON_DIR)+config_path);
    mDfobj->Initialize(mSoftWorld);

    printf("-------- 1\n");
    // mDfcloth = new Dfcloth<TinyScalar, TinyConstants>(5E5,2E5,0.4);
    // mDfcloth->SetDfobj(std::string(SOFTCON_DIR)+"/data/json_files/sphere_3_layers.json");
    // mDfcloth->Initialize(mSoftWorld);

    

    printf("-------- 2\n");
    mSoftWorld->Initialize();
    printf("-------- 3\n");
    mDfobj->SetInitReferenceRotation(mSoftWorld->GetPositions());

    printf("-------- 4\n");
    InitializeActions();

    printf("-------- 5\n");
    // mTargetVelocity.setZero();

    // mAverageVelocity.setZero();
    // mAverageVelocityDeque.clear();

    // printf("-------- 6\n");
    // UpdateRandomTargetVelocity();
    printf("-------- 7\n");
}
// template <typename TinyScalar, typename TinyConstants> 
// void Environment<TinyScalar, TinyConstants>::
// Step()
// {
//     static int total_step = 0;
//     std::cout << "========================================Time Step: " << total_step << "\n";
//     if (total_step==54)
//         exit(0);
//     // Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> external_force = mDfobj->ComputeDragForces(mSoftWorld);
//     Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> external_force = mSoftWorld->mExternalForces;
//     external_force.setZero();
//     // for (int i = 0; i < external_force.size() / 3; i++) {
//     // for (int i = 0; i < 100; i++) {
//     //         external_force[3*i+2] = -1;
//     // }
//     mSoftWorld->SetExternalForce(external_force);

//     // for (int i = 0; i < 10; i++)
//     //     printf("![%f %f] ", external_force[i*3+2], mSoftWorld->mExternalForces[i*3+2]);
//     mSoftWorld->TimeStepping();

//     // Eigen::Matrix<TinyScalar, 3, 1> v_front = Eigen::Matrix<TinyScalar, 3, 1>(0,0,0);
//     // if(mAverageVelocityDeque.size() > mSimulationHz) {
//     //     v_front = mAverageVelocityDeque.front();
//     //     mAverageVelocityDeque.pop_front();
//     // }
//     // printf("---------- 5\n");
//     // const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& tmp_vs = mSoftWorld->GetVelocities();
//     // Eigen::Matrix<TinyScalar, 3, 1> v_center = tmp_vs.template block<3,1>(3*mDfobj->GetCenterIndex(),0);
//     // // Eigen::Matrix<TinyScalar, 3, 1> v_center ;
//     // // mSoftWorld->GetVelocities().template block<3,1>(3*mDfobj->GetCenterIndex(),0);
//     // mAverageVelocityDeque.push_back(v_center);
//     // mAverageVelocity = mAverageVelocity - (v_front)/mSimulationHz + v_center/mSimulationHz;

//     char filename[50];
//     sprintf(filename, "out/%03d.obj",total_step++);
//     OutputSurface(filename);
// }


// template <typename TinyScalar, typename TinyConstants> 
// void Environment<TinyScalar, TinyConstants>::
// Step()
// {
//     // Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> external_force = mDfobj->ComputeDragForces(mSoftWorld);
//     // mSoftWorld->SetExternalForce(external_force);

//     mSoftWorld->TimeSteppingRigid();

//     Eigen::Matrix<TinyScalar, 3, 1> v_front = Eigen::Matrix<TinyScalar, 3, 1>(0,0,0);
//     if(mAverageVelocityDeque.size() > mSimulationHz) {
//         v_front = mAverageVelocityDeque.front();
//         mAverageVelocityDeque.pop_front();
//     }
//     const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& tmp_vs = mSoftWorld->GetVelocities();
//     Eigen::Matrix<TinyScalar, 3, 1> v_center = tmp_vs.template block<3,1>(3*mDfobj->GetCenterIndex(),0);
//     // Eigen::Matrix<TinyScalar, 3, 1> v_center ;
//     // mSoftWorld->GetVelocities().template block<3,1>(3*mOctopus->GetCenterIndex(),0);
//     mAverageVelocityDeque.push_back(v_center);
//     mAverageVelocity = mAverageVelocity - (v_front)/mSimulationHz + v_center/mSimulationHz;
// }

template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
Step()
{
    // static int total_step = 0;
    // std::cout << "========================================Time Step: " << total_step++ << "\n";
    // Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> act = mSoftWorld->mJointTorque;
    // act.fill(TinyConstants::scalar_from_double(0.));
    // std::cout << "action size " << act.rows() << std::endl;
    // for (int i = 0; i < act.rows(); ++i)
    //     act[i]= 6* sin( 2./30. * total_step * 3.1415926535 );
    // SetActions(act);
    
    // if (total_step==54)
    //     exit(0);
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> external_force = mDfobj->ComputeDragForces(mSoftWorld);
    mSoftWorld->SetExternalForce(external_force);
    // external_force.setZero();

    mSoftWorld->TimeStepping();

    Eigen::Matrix<TinyScalar, 3, 1> v_front = Eigen::Matrix<TinyScalar, 3, 1>(0,0,0);
    if(mAverageVelocityDeque.size() > mSimulationHz) {
        v_front = mAverageVelocityDeque.front();
        mAverageVelocityDeque.pop_front();
    }
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tmp_vs = mSoftWorld->GetVelocities();
    Eigen::Matrix<TinyScalar, 3, 1> v_center;
    v_center = tmp_vs.template block<3,1>(3*mDfobj->GetCenterIndex(),0);
    // Eigen::Matrix<TinyScalar, 3, 1> v_center ;
    // mSoftWorld->GetVelocities().template block<3,1>(3*mOctopus->GetCenterIndex(),0);
    mAverageVelocityDeque.push_back(v_center);
    mAverageVelocity = mAverageVelocity - (v_front)/mSimulationHz + v_center/mSimulationHz;


    // const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x1 = GetSoftWorld()->GetPositions();
    // char filename[50];
    // sprintf(filename, "out/%03d.obj",total_step++);
    // GetDfobj()->OutputSurface(filename, x1);

}


template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
Step(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x,
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& v,
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& jq,
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& jqd,
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& tau,
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& next_x,
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& next_v,
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& next_jq,
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& next_jqd)
{
    // static int total_step = 0;
    // std::cout << "========================================Time Step: " << total_step++ << "\n";
    mSoftWorld->mX = x;
    mSoftWorld->mV = v;
    mSoftWorld->mJointQ = jq;
    mSoftWorld->mJointQd = jqd;



    // external_force.setZero();
    if (!mDfobj->mIsUnderwater) {
        Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> g_total = 
            Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>::Zero(3*mSoftWorld->GetNumVertices());
        for (int i = 0; i < mSoftWorld->GetNumVertices(); i++)
            g_total[i*3+1] += mDfobj->mGravityY;
        mSoftWorld->mV += g_total*mSoftWorld->mTimeStep;
    } else {
        Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> external_force = mDfobj->ComputeDragForces(mSoftWorld);
        mSoftWorld->SetExternalForce(external_force);
    }

    mSoftWorld->TimeStepping();
    next_x = mSoftWorld->mX;
    next_v = mSoftWorld->mV;
    next_jq = mSoftWorld->mJointQ;
    next_jqd = mSoftWorld->mJointQd;


    // const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x1 = next_x;
    // char filename[50];
    // sprintf(filename, "out/%03d.obj",total_step++);
    // GetDfobj()->OutputSurface(filename, x1);
}

template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
SaveObj(std::string filename, const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &x) {
    const auto& soft_solid = GetDfobj();
    FILE *file = fopen(filename.c_str(), "w");
    const auto& vertices = soft_solid->mMesh->GetVertices();
    for (int i = 0; i < vertices.size(); ++i)
        fprintf(file, "v %06f %06f %06f\n",
            TinyConstants::getDouble(x[i*3+0]),
            TinyConstants::getDouble(x[i*3+1]),
            TinyConstants::getDouble(x[i*3+2]));
    for (int i = 0; i < soft_solid->mContours.size(); ++i)
        fprintf(file, "f %d %d %d\n",
            soft_solid->mContours[i][0]+1,
            soft_solid->mContours[i][1]+1,
            soft_solid->mContours[i][2]+1);
    fclose(file);
}

// template <typename TinyScalar, typename TinyConstants> 
// void Environment<TinyScalar, TinyConstants>::
// Step()
// {
//     printf("---------- 2\n");
//     // Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> external_force = mDfobj->ComputeDragForces(mSoftWorld);
//     // mSoftWorld->SetExternalForce(external_force);

//     printf("---------- 3\n");
//     mSoftWorld->TimeStepping();

//     printf("---------- 4\n");
//     // Eigen::Matrix<TinyScalar, 3, 1> v_front = Eigen::Matrix<TinyScalar, 3, 1>(0,0,0);
//     // if(mAverageVelocityDeque.size() > mSimulationHz) {
//     //     v_front = mAverageVelocityDeque.front();
//     //     mAverageVelocityDeque.pop_front();
//     // }
//     // printf("---------- 5\n");
//     // const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& tmp_vs = mSoftWorld->GetVelocities();
//     // Eigen::Matrix<TinyScalar, 3, 1> v_center = tmp_vs.template block<3,1>(3*mDfobj->GetCenterIndex(),0);
//     // // Eigen::Matrix<TinyScalar, 3, 1> v_center ;
//     // // mSoftWorld->GetVelocities().template block<3,1>(3*mDfobj->GetCenterIndex(),0);
//     // mAverageVelocityDeque.push_back(v_center);
//     // mAverageVelocity = mAverageVelocity - (v_front)/mSimulationHz + v_center/mSimulationHz;
// }

template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
OutputSurface(std::string filename) {
    FILE *file = fopen(filename.c_str(), "w");
    const auto& vertices = mSoftWorld->mCollision_X_next;
    const auto& triangles = mSoftWorld->mCollisionMesh->m_triangles;

    for (int i = 0; i < vertices.size()/3; ++i)
        fprintf(file, "v %06f %06f %06f\n",
            vertices[i*3+0],vertices[i*3+1],vertices[i*3+2]);
    for (int i = 0; i < triangles.size()/3; ++i)
        fprintf(file, "f %d %d %d\n",
            triangles[i*3+0]+1,triangles[i*3+1]+1,triangles[i*3+2]+1);
    fclose(file);
}

template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
SetPhase(const int& phase)
{
    mPhase = phase;
}

template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
Reset()
{
    mAverageVelocityDeque.clear();
    mAverageVelocity.setZero();
    for(int i=0;i<mDfobj->GetMuscles().size();i++)
    {
        mDfobj->GetMuscles()[i]->Reset();
    }
    mSoftWorld->Reset();
    mPhase = 0;
}
// DeepRL

template <typename TinyScalar, typename TinyConstants> 
const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& Environment<TinyScalar, TinyConstants>::
GetStates()
{
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x = mSoftWorld->GetPositions();
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& v = mSoftWorld->GetVelocities();
    Eigen::Matrix<TinyScalar, 3, 1> center_position;
    center_position = x.template block<3,1>(3*mDfobj->GetCenterIndex(),0);
    
    Eigen::Matrix<TinyScalar, 3, 3> R = mDfobj->GetReferenceRotation(LOCAL,x);

    Eigen::Matrix<TinyScalar, 3, 1> local_target_velocity = R*mTargetVelocity;
    Eigen::Matrix<TinyScalar, 3, 1> local_average_velocity = R*mAverageVelocity;

    const std::vector<int>& sampling_index = mDfobj->GetSamplingIndex();
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> local_position(3*sampling_index.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> local_velocity(3*sampling_index.size());

    for(int i=0; i<sampling_index.size(); i++) {
        local_position.template block<3,1>(3*i,0) = R*x.template block<3,1>(3*sampling_index[i],0);
        local_velocity.template block<3,1>(3*i,0) = R*v.template block<3,1>(3*sampling_index[i],0);
    }

    int num_states = local_target_velocity.size() + local_average_velocity.size()
        + local_position.size() + local_velocity.size();

    mStates.resize(num_states);
    mStates.setZero();

    mStates<<local_target_velocity,local_average_velocity,
        local_position,local_velocity;

    return mStates;
}

template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
InitializeActions()
{
    const auto& muscles = mDfobj->GetMuscles();

    int num_action =4*muscles.size();
    mActions.resize(num_action);
    mActions.setZero();

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> real_lower_bound(num_action);
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> real_upper_bound(num_action);

    std::cout << num_action << " num_action\n";
    std::cout << muscles.size() << " muscles\n";
    if (muscles.size() == 0)
        return;
    int cnt =0;
    for(const auto& m : muscles) 
    {
        real_lower_bound.segment(cnt,4) = m->GetActionLowerBound();
        real_upper_bound.segment(cnt,4) = m->GetActionUpperBound();
        cnt+=4;
    }

    mNormLowerBound.resize(real_lower_bound.size());
    mNormUpperBound.resize(real_upper_bound.size());

    mNormLowerBound.setOnes();
    mNormLowerBound *= -5.0;
    mNormUpperBound.setOnes();
    mNormUpperBound *= 5.0;

    mNormalizer = new Normalizer<TinyScalar, TinyConstants>(real_upper_bound,real_lower_bound,
        mNormUpperBound,mNormLowerBound);
}

template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
SetActions(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& actions)
{
    if (mDfobj->mIsSkeletonAct) {
        mSoftWorld->mJointTorque = actions;
    } else {
        Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> real_actions = mNormalizer->NormToReal(actions);
        mDfobj->SetActivationLevels(real_actions,mPhase);
    }
    // mPhase+=1;

}


template <typename TinyScalar, typename TinyConstants> 
std::map<std::string,TinyScalar> Environment<TinyScalar, TinyConstants>::
GetRewards()
{
    std::map<std::string,TinyScalar> reward_map;

    TinyScalar d = (mAverageVelocity-mTargetVelocity).norm();
    TinyScalar reward_target = exp(-d*d/0.05);

    Eigen::Matrix<TinyScalar, 3, 1> v_face = mDfobj->GetForwardVector(mSoftWorld->GetPositions());
    Eigen::Matrix<TinyScalar, 3, 1> v_tar_dir = mTargetVelocity.normalized();
    auto d_direction = (v_face - v_tar_dir).norm();
    TinyScalar reward_direction = exp(-fabs(d_direction)/0.1);

    TinyScalar w_target = 1.0;
    TinyScalar w_direction = 2.0;

    TinyScalar reward = w_target*reward_target + w_direction*reward_direction;

    reward_map["target"] = w_target*reward_target;
    reward_map["direction"] = w_direction*reward_direction;
    reward_map["total"] = reward;

    return reward_map;
}

template <typename TinyScalar, typename TinyConstants> 
bool Environment<TinyScalar, TinyConstants>::
isEndOfEpisode()
{
    bool eoe =false; 

    return eoe;
}

template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
UpdateRandomTargetVelocity()
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> l_distribution(0.75,0.75);
	std::uniform_real_distribution<double> angle_distribution(0.0,0.0);
	std::uniform_real_distribution<double> vx_distribution(-1.0,1.0);
	std::uniform_real_distribution<double> vy_distribution(-1.0,1.0);
	std::uniform_real_distribution<double> vz_distribution(-1.0,1.0);


    Eigen::Matrix<TinyScalar, 3, 1> v,axis;

    TinyScalar length = TinyConstants::scalar_from_double(l_distribution(generator));

    v = length*mDfobj->GetForwardVector(mSoftWorld->GetPositions()).normalized();

    // TinyScalar angle = angle_distribution(generator);
    // axis[0] = vx_distribution(generator);
    // axis[1] = vy_distribution(generator);
    // axis[2] = vz_distribution(generator);

    // Eigen::Matrix<TinyScalar, 3, 3> R;
    // R = Eigen::AngleAxisd(angle,axis.normalized());

    // mTargetVelocity = R*v;
    mTargetVelocity = v;
}
#endif