#ifndef __ENVIRONMENT_H__
#define __ENVIRONMENT_H__
#include <deque>
#include "fem/World.h"
#include "Dfobj.h"
#include "Dfcloth.h"
#include "Dfobstacle.h"
#include "Normalizer.h"
// class Dfobj;
// class Normalizer;

template <typename TinyScalar, typename TinyConstants> 
class Environment
{
public:
    Environment(int total_step, const std::string& config_path);

	FEM::World<TinyScalar, TinyConstants>* GetSoftWorld(){return mSoftWorld;};
    std::vector<Dfobj<TinyScalar, TinyConstants>*>& GetDfobjs(){return mDfobjs;};
    std::vector<Dfobstacle<TinyScalar, TinyConstants>*>& GetDfobss(){return mDfobstacles;};
    std::vector<Dfcloth<TinyScalar, TinyConstants>*>& GetDfcloths(){return mDfcloths;};

	TinyScalar GetSimulationHz(){return mSimulationHz;};

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
    void Save(std::string filename, const Dfobj<TinyScalar, TinyConstants> *soft_solid, const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &x);
    void SaveObj(std::string filename, const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &x);

	void Reset();

	// const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& GetStates();

	// void InitializeActions();
	// const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& GetActions() {return mActions;};
	void SetActions(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& actions);

	std::map<std::string,TinyScalar> GetRewards();

    bool isEndOfEpisode(){return false;};

	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> GetNormLowerBound() {assert(false);};
	Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> GetNormUpperBound() {assert(false);};
	Normalizer<TinyScalar, TinyConstants>*	GetNormalizer() {assert(false);};

    void IncPhase();
	void SetPhase(const int& phase); 
	const int& GetPhase() {return mPhase;};

	// Eigen::Matrix<TinyScalar, 3, 1> GetAverageVelocity() {return mAverageVelocity;};
	Eigen::Matrix<TinyScalar, 3, 1> GetTargetVelocity() {assert(false);};
	// void UpdateRandomTargetVelocity();
    void OutputSurface(std::string filename);

	FEM::World<TinyScalar, TinyConstants>*					mSoftWorld;
    // Dfcloth<TinyScalar, TinyConstants>*                       mDfcloth;
    // Dfobj<TinyScalar, TinyConstants>*                       mDfobj;
    std::vector<Dfcloth<TinyScalar, TinyConstants>*>           mDfcloths;
    std::vector<Dfobstacle<TinyScalar, TinyConstants>*>           mDfobstacles;
    std::vector<Dfobj<TinyScalar, TinyConstants>*>             mDfobjs;

	int 							mSimulationHz;
    int                             mTotalStep;

	int 							mPhase;
    bool mIsUnderwater;
    Eigen::Matrix<TinyScalar, 3, 1> mGravity;
};




#include <chrono>
#include <random>
#include <fstream>

template <typename TinyScalar, typename TinyConstants> 
Environment<TinyScalar, TinyConstants>::
Environment(int total_step, const std::string& config_path)
    :mTotalStep(total_step)
{
    string path = std::string(SOFTCON_DIR)+config_path;
    std::cout <<"+ SetDfobj Data Path: " << path << std::endl;
    std::FILE* file = std::fopen(path.c_str(), "r");
    if (file == nullptr)
        throw std::invalid_argument(std::string("Couldn't open file ") + path);
    constexpr std::size_t buffer_size = 10240;
    char buffer[buffer_size];

    printf("loading scene\n");
    rapidjson::FileReadStream is = rapidjson::FileReadStream(file, buffer, buffer_size);
    rapidjson::Document scene_description;
    scene_description.ParseStream(is);
    fclose(file);
    printf("finish scene\n");

    if (!scene_description.IsObject())
        throw std::invalid_argument("JSON reader - Error : root must be an object");

    printf("loading Parameters\n");
    const rapidjson::Value& para_array = scene_description["Parameters"];
    const auto& para_object = para_array[0];
    if (para_object.HasMember("Simulation HZ"))
        mSimulationHz = para_object["Simulation HZ"].GetInt();
    else
        mSimulationHz = 200;

    bool do_collision = false;
    if (para_object.HasMember("Handle collision"))
        do_collision = para_object["Handle collision"].GetBool();

    int max_iteration = 30;
    if (para_object.HasMember("Max Iteration"))
        max_iteration = para_object["Max Iteration"].GetInt();

    TinyScalar damping_coeff = 0.9999;
    if (para_object.HasMember("Damping coefficient"))
        damping_coeff = para_object["Damping coefficient"].GetDouble();



    std::cout << mTotalStep << " mTotalStep\n";
    std::cout << config_path << " config_path\n";
    std::cout << mSimulationHz << " simulationHz\n";
    std::cout << do_collision << " Handle collision\n";


    mSoftWorld = new FEM::World<TinyScalar, TinyConstants>(
        1.0/mSimulationHz,
        max_iteration,
        damping_coeff,
        do_collision
    );
    if (para_object.HasMember("World Alpha"))
        mSoftWorld->mAlpha = para_object["World Alpha"].GetDouble();
    else
        mSoftWorld->mAlpha = 1;
    if (para_object.HasMember("World Friction"))
        mSoftWorld->mFriction = para_object["World Friction"].GetDouble();
    else
        mSoftWorld->mFriction = 1;
    if (para_object.HasMember("Max Iteration"))
        mSoftWorld->mMaxIteration1 = para_object["Max Iteration"].GetInt();
    else
        mSoftWorld->mMaxIteration1 = 30;
    if (para_object.HasMember("Self-collision"))
        mSoftWorld->m_handle_self_collision = para_object["Self-collision"].GetBool();


    mDfobjs.clear();
    mDfcloths.clear();
    mDfobstacles.clear();

    if (para_object.HasMember("Is underwater")) {
        printf("loading Is underwater\n");
        mIsUnderwater = para_object["Is underwater"].GetInt();
    } else {
        mIsUnderwater = 0;
    }

    if (para_object.HasMember("Gravity")) {
        printf("loading Gravity\n");

        const auto &gravity = para_object["Gravity"];
        if (gravity.Size() != 3)
            std::cout << "The Gravity vector should have the length of 3 but got length " 
                        << gravity.Size() << "\n";
        for (int i = 0; i < 3; ++i)
            mGravity[i] = gravity[i].GetDouble();
    } else {
        for (int i = 0; i < 3; ++i)
            mGravity[i] = 0.;
    }
    
    int numdfobjs = scene_description["Meshes"].Size();
    for (int i = 0; i < numdfobjs; ++i) {
        auto mDfobj = new Dfobj<TinyScalar, TinyConstants>(5E5,2E5,0.4);
        mDfobj->SetDfobj(
            scene_description["Meshes"][i], 
            Eigen::Matrix<TinyScalar, 3, 1>(0.,0.,0.));
        mDfobj->Initialize(mSoftWorld);
        mDfobjs.push_back(mDfobj);
    }

    numdfobjs = scene_description["Clothes"].Size();
    for (int i = 0; i < numdfobjs; ++i) {
        auto mDfcloth = new Dfcloth<TinyScalar, TinyConstants>(5E5,2E5,0.4);
        mDfcloth->SetDfobj(
            scene_description["Clothes"][i]);
        mDfcloth->Initialize(mSoftWorld);
        mDfcloths.push_back(mDfcloth);
    }

    numdfobjs = scene_description["Obstacles"].Size();
    for (int i = 0; i < numdfobjs; ++i) {
        auto mDfobstacle = new Dfobstacle<TinyScalar, TinyConstants>();
        mDfobstacle->SetDfobj(
            scene_description["Obstacles"][i], 
            Eigen::Matrix<TinyScalar, 3, 1>(0.,0.,0.));
        mDfobstacle->Initialize(mSoftWorld);
        mDfobstacles.push_back(mDfobstacle);
    }

    printf("-------- 1\n");

    mSoftWorld->Initialize();
    int offset = 0;
    int num_action = 0;
    for (auto mDfobj : mDfobjs) {
        mDfobj->SetInitReferenceRotation(mSoftWorld->GetPositions().segment(offset, mDfobj->mNumUnknown));
        mDfobj->InitializeActions();
        offset += mDfobj->mNumUnknown;
        num_action += mDfobj->mActions.size();
    }
    // mActions.resize(num_action);
    mPhase = 0;
}


template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
Step()
{   
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> external_force = mSoftWorld->GetPositions();
    external_force.setZero();
    int offset = 0;
    for (auto mDfobj : mDfobjs) {
        if (!mIsUnderwater) {
            if (mDfobj->mEnableGravity) {
                Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> g_total = 
                    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>::Zero(mDfobj->mNumUnknown);
                for (int i = 0; i < mDfobj->mNumUnknown/3; i++) 
                    g_total.segment(i * 3, 3) += mGravity;
                    // for (int j = 0; j < 3; j++)
                    //     g_total[i * 3 + j] += mGravity[j];
                mSoftWorld->mV.segment(offset, mDfobj->mNumUnknown) += g_total*mSoftWorld->mTimeStep;
            } 
        } else {
            external_force.segment(offset, mDfobj->mNumUnknown) = mDfobj->ComputeDragForces(mSoftWorld);
            
            // std::cout <<"external_force "<<external_force[0]<<" "<<external_force[1]<<" "<<external_force[2]<<" "<<external_force[3]<<std::endl;
        }
        offset += mDfobj->mNumUnknown;
    }
    mSoftWorld->SetExternalForce(external_force);

    mSoftWorld->TimeStepping();



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

    Step();

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
Save(std::string filename, const Dfobj<TinyScalar, TinyConstants> *soft_solid, const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &x) {
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

template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
SaveObj(std::string filename, const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &x) {
    int offset = 0;
    for (int i = 0; i < mDfobjs.size(); ++i) {
        std::string save_file_name = filename + "_DfObj" + std::to_string(i) + ".obj";
        Save(save_file_name, mDfobjs[i], x.segment(offset, mDfobjs[i]->mNumUnknown));
        offset += mDfobjs[i]->mNumUnknown;
    }
    for (int i = 0; i < mDfcloths.size(); ++i) {
        std::string save_file_name = filename + "_DfCloths" + std::to_string(i) + ".obj";
        Save(save_file_name, mDfcloths[i], x.segment(offset, mDfcloths[i]->mNumUnknown));
        offset += mDfcloths[i]->mNumUnknown;
    }
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
IncPhase()
{
    mPhase += 1;
    for (auto mDfobj : mDfobjs)
        mDfobj->mPhase = mPhase;
}

template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
SetPhase(const int& phase)
{
    for (auto mDfobj : mDfobjs)
        mDfobj->mPhase = phase;
}

template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
Reset()
{
    for (auto mDfobj : mDfobjs) {
        for(int i=0;i<mDfobj->GetMuscles().size();i++)
        {
            mDfobj->GetMuscles()[i]->Reset();
        }
        mDfobj->mPhase = 0;
    }
    mSoftWorld->Reset();
}
// DeepRL

// template <typename TinyScalar, typename TinyConstants> 
// const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& Environment<TinyScalar, TinyConstants>::
// GetStates()
// {
//     const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x = mSoftWorld->GetPositions();
//     const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& v = mSoftWorld->GetVelocities();
//     Eigen::Matrix<TinyScalar, 3, 1> center_position;
//     center_position = x.template block<3,1>(3*mDfobj->GetCenterIndex(),0);
    
//     Eigen::Matrix<TinyScalar, 3, 3> R = mDfobj->GetReferenceRotation(LOCAL,x);

//     Eigen::Matrix<TinyScalar, 3, 1> local_target_velocity = R*mTargetVelocity;
//     Eigen::Matrix<TinyScalar, 3, 1> local_average_velocity = R*mAverageVelocity;

//     const std::vector<int>& sampling_index = mDfobj->GetSamplingIndex();
//     Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> local_position(3*sampling_index.size());
//     Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> local_velocity(3*sampling_index.size());

//     for(int i=0; i<sampling_index.size(); i++) {
//         local_position.template block<3,1>(3*i,0) = R*x.template block<3,1>(3*sampling_index[i],0);
//         local_velocity.template block<3,1>(3*i,0) = R*v.template block<3,1>(3*sampling_index[i],0);
//     }

//     int num_states = local_target_velocity.size() + local_average_velocity.size()
//         + local_position.size() + local_velocity.size();

//     mStates.resize(num_states);
//     mStates.setZero();

//     mStates<<local_target_velocity,local_average_velocity,
//         local_position,local_velocity;

//     return mStates;
// }

// template <typename TinyScalar, typename TinyConstants> 
// void Environment<TinyScalar, TinyConstants>::
// InitializeActions()
// {
//     const auto& muscles = mDfobj->GetMuscles();

//     int num_action =4*muscles.size();
//     mActions.resize(num_action);
//     mActions.setZero();

//     Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> real_lower_bound(num_action);
//     Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> real_upper_bound(num_action);

//     std::cout << num_action << " num_action\n";
//     std::cout << muscles.size() << " muscles\n";
//     if (muscles.size() == 0)
//         return;
//     int cnt =0;
//     for(const auto& m : muscles) 
//     {
//         real_lower_bound.segment(cnt,4) = m->GetActionLowerBound();
//         real_upper_bound.segment(cnt,4) = m->GetActionUpperBound();
//         cnt+=4;
//     }

//     mNormLowerBound.resize(real_lower_bound.size());
//     mNormUpperBound.resize(real_upper_bound.size());

//     mNormLowerBound.setOnes();
//     mNormLowerBound *= -5.0;
//     mNormUpperBound.setOnes();
//     mNormUpperBound *= 5.0;

//     mNormalizer = new Normalizer<TinyScalar, TinyConstants>(real_upper_bound,real_lower_bound,
//         mNormUpperBound,mNormLowerBound);
// }

template <typename TinyScalar, typename TinyConstants> 
void Environment<TinyScalar, TinyConstants>::
SetActions(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& actions)
{
    int act_now = 0;
    for (auto mDfobj : mDfobjs) {
        mDfobj->SetActions(
            actions.segment(act_now, mDfobj->mActions.size()));
        act_now += mDfobj->mActions.size();
    }

}


// template <typename TinyScalar, typename TinyConstants> 
// std::map<std::string,TinyScalar> Environment<TinyScalar, TinyConstants>::
// GetRewards()
// {
//     std::map<std::string,TinyScalar> reward_map;

//     TinyScalar d = (mAverageVelocity-mTargetVelocity).norm();
//     TinyScalar reward_target = exp(-d*d/0.05);

//     Eigen::Matrix<TinyScalar, 3, 1> v_face = mDfobj->GetForwardVector(mSoftWorld->GetPositions());
//     Eigen::Matrix<TinyScalar, 3, 1> v_tar_dir = mTargetVelocity.normalized();
//     auto d_direction = (v_face - v_tar_dir).norm();
//     TinyScalar reward_direction = exp(-fabs(d_direction)/0.1);

//     TinyScalar w_target = 1.0;
//     TinyScalar w_direction = 2.0;

//     TinyScalar reward = w_target*reward_target + w_direction*reward_direction;

//     reward_map["target"] = w_target*reward_target;
//     reward_map["direction"] = w_direction*reward_direction;
//     reward_map["total"] = reward;

//     return reward_map;
// }

// template <typename TinyScalar, typename TinyConstants> 
// bool Environment<TinyScalar, TinyConstants>::
// isEndOfEpisode()
// {
//     bool eoe =false; 

//     return eoe;
// }

// template <typename TinyScalar, typename TinyConstants> 
// void Environment<TinyScalar, TinyConstants>::
// UpdateRandomTargetVelocity()
// {
// 	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
// 	std::default_random_engine generator(seed);
// 	std::uniform_real_distribution<double> l_distribution(0.75,0.75);
// 	std::uniform_real_distribution<double> angle_distribution(0.0,0.0);
// 	std::uniform_real_distribution<double> vx_distribution(-1.0,1.0);
// 	std::uniform_real_distribution<double> vy_distribution(-1.0,1.0);
// 	std::uniform_real_distribution<double> vz_distribution(-1.0,1.0);


//     Eigen::Matrix<TinyScalar, 3, 1> v,axis;

//     TinyScalar length = TinyConstants::scalar_from_double(l_distribution(generator));

//     v = length*mDfobj->GetForwardVector(mSoftWorld->GetPositions()).normalized();

//     // TinyScalar angle = angle_distribution(generator);
//     // axis[0] = vx_distribution(generator);
//     // axis[1] = vy_distribution(generator);
//     // axis[2] = vz_distribution(generator);

//     // Eigen::Matrix<TinyScalar, 3, 3> R;
//     // R = Eigen::AngleAxisd(angle,axis.normalized());

//     // mTargetVelocity = R*v;
//     mTargetVelocity = v;
// }

#endif