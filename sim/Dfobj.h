#ifndef __Dfskeleton_H__
#define __Dfskeleton_H__
#include "fem/World.h"
#include "fem/Mesh/MeshHeaders.h"
#include "Normalizer.h"
#include "Muscle.h"

#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#include <fstream>
#include <math.h>
#define LOCAL   0
#define GLOBAL  1

template <typename TinyScalar, typename TinyConstants> 
class Dfobj
{
public:
    Dfobj(const TinyScalar& muscle_stiffness,const TinyScalar& youngs_modulus,const TinyScalar& poisson_ratio);
    ~Dfobj();

    virtual void Initialize(FEM::World<TinyScalar, TinyConstants>* world);
    virtual void SetDfobj(const rapidjson::Value &object, Eigen::Matrix<TinyScalar, 3, 1> offset);
    void SetMesh(const std::string& path,const Eigen::Transform<TinyScalar, 3, 2>& T);
    void MakeMuscles(const std::string& path,const TinyScalar& gamma);

    void SetActivationLevels(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& actions,const int& phase); 
    void SetKey(std::string filename);

    Eigen::Matrix<TinyScalar, 3, 3> GetReferenceRotation(bool type,const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
    void SetInitReferenceRotation(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);

    const std::vector<Eigen::Vector3i>& GetContours(){return mContours;};

    virtual void SetVertexNormal();
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& GetVertexNormal(){return mVertexNormal;};

    const std::vector<Muscle<TinyScalar, TinyConstants>*>& GetMuscles() {return mMuscles;};

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ComputeDragForces(FEM::World<TinyScalar, TinyConstants>* world);

    Eigen::Matrix<TinyScalar, 3, 1> GetUpVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
    Eigen::Matrix<TinyScalar, 3, 1> GetForwardVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x); 

    const int& GetCenterIndex() {return mCenterIndex;};
    const int& GetEndEffectorIndex() {return mEndEffectorIndex;};
    const std::vector<int>& GetSamplingIndex() {return mSamplingIndex;};  
    void OutputSurface(std::string filename, const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &x);
    void InitializeActions();
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& GetActions() {return mActions;};
    void SetActions(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& actions);

    void SetActivationLevelsAggregate(Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& v);
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> GetActivationLevelsAggregate();
// private:
    TinyScalar                                  mMuscleStiffness;
    TinyScalar                                  mYoungsModulus;
    TinyScalar                                  mPoissonRatio;

    FEM::World<TinyScalar, TinyConstants>*      mSoftWorld;
    FEM::Mesh<TinyScalar, TinyConstants>*       mMesh;
    Eigen::Transform<TinyScalar, 3, 2>                         mScalingMatrix;

    std::vector<FEM::Constraint<TinyScalar, TinyConstants>*>            mConstraints;

    std::vector<Muscle<TinyScalar, TinyConstants>*>                 mMuscles;

    std::vector<Eigen::Vector3i>            mContours;

    int                                     mCenterIndex;
    int                                     mNumUnknown = 0;
    int                                     mEndEffectorIndex;
    bool                                    mEnableGravity;

    std::vector<int>                        mSamplingIndex;
    Eigen::Vector3i                         mLocalContourIndex;

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                         mVertexNormal;

    Eigen::Matrix<TinyScalar, 3, 3>                         mInitReferenceMatrix;

    bool                                    mIsVolumeAct;
    int                                     mIsSkeletonAct, mActuatorType;
    TinyScalar                              mAttachStiffness = 0.;
    bool                                    mIsAttach = false;
    // std::vector<FEM::CorotateFEMConstraintActuate<TinyScalar, TinyConstants>*> mVolumeMuscles;
    std::vector<FEM::PneumaticConstraint<TinyScalar, TinyConstants>*> mVolumeMuscles;


    std::vector<int>                    mContactIdx;
    std::vector<Eigen::Vector3i>        mContactFace;
    std::vector<TinyScalar> m_masses;

    TinyScalar C_d(TinyScalar theta);
    TinyScalar C_t(TinyScalar theta);
    TinyScalar Cd(TinyScalar aoa,TinyScalar a,TinyScalar c); 
    TinyScalar _interpolate(TinyScalar t, TinyScalar a, TinyScalar b);
    TinyScalar Cl(TinyScalar aoa, TinyScalar cutoff, TinyScalar x[5], TinyScalar y[5]); 

    std::vector<int>                        mnonrigid, mrigid;
    std::vector<std::vector<int>>           mRigidBodies;
    std::vector<Link<TinyScalar, TinyConstants>*>                      mlinks;

    Normalizer<TinyScalar, TinyConstants>*                  mNormalizer;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                 mNormLowerBound;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                 mNormUpperBound;

    // Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                  mStates;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                 mActions;

    int                             mPhase;
    std::set<int>                   mVolumeVertIdx, mAttachIdx;
    TinyScalar                      mTotalMass = 10.0;
    bool                            mIsMassSame = true;
};




using namespace FEM;

template <typename TinyScalar, typename TinyConstants> 
Dfobj<TinyScalar, TinyConstants>::
Dfobj(const TinyScalar& muscle_stiffness,const TinyScalar& youngs_modulus,const TinyScalar& poisson_ratio)
    :mMesh()
    ,mMuscleStiffness(muscle_stiffness)
    ,mYoungsModulus(youngs_modulus)
    ,mPhase(0)
    ,mPoissonRatio(poisson_ratio),mIsVolumeAct(true)
{
    printf("construct Dfskeleton\n");
}


template <typename TinyScalar, typename TinyConstants> 
Dfobj<TinyScalar, TinyConstants>::
~Dfobj() {
    delete mMesh;
}

template <typename TinyScalar, typename TinyConstants> 
void Dfobj<TinyScalar, TinyConstants>::
SetDfobj(const rapidjson::Value &object, Eigen::Matrix<TinyScalar, 3, 1> offset)
{
    if (object.HasMember("Is skeleton actuated")) {
        printf("loading Is skeleton actuated\n");
        mIsSkeletonAct = object["Is skeleton actuated"].GetInt();
    } else {
        mIsSkeletonAct = 0;
    }

    // 0 no actuator
    // 1 skeleton
    // 2 muscle
    // 3 pneumatic
    if (object.HasMember("Actuator type")) {
        printf("loading Actuator type\n");
        mActuatorType = object["Actuator type"].GetInt();
    } else {
        mActuatorType = 0;
    }

    assert(object.HasMember("Obj filename"));
    printf("loading mesh\n");
    double scale = 0.05;
    if (object.HasMember("Scale"))
        scale = object["Scale"].GetDouble();
    std::cout << scale << " scale\n" ;
    assert(scale > 0);
    if (object.HasMember("Enable Gravity"))
        mEnableGravity = object["Enable Gravity"].GetInt();
    else
        mEnableGravity = true;
    if (object.HasMember("Youngs Modulus"))
        mYoungsModulus = object["Youngs Modulus"].GetDouble();
    if (object.HasMember("Muscle Stiffness"))
        mMuscleStiffness = object["Muscle Stiffness"].GetDouble();
    if (object.HasMember("Poisson Ratio"))
        mPoissonRatio = object["Poisson Ratio"].GetDouble();
    if (object.HasMember("Mass"))
        mTotalMass = object["Mass"].GetDouble();
    if (object.HasMember("Is Mass Same"))
        mIsMassSame = object["Is Mass Same"].GetBool();

    

    std::string rest_state_file_name = object["Obj filename"].GetString();
    Eigen::Transform<TinyScalar, 3, 2> T=Eigen::Transform<TinyScalar, 3, 2>::Identity();
    T(0,0) *= scale; T(1,1) *= scale; T(2,2) *= scale; 
    T.translation() = offset;
    SetMesh(rest_state_file_name,T);



    if (object.HasMember("Muscle filename")) {
        printf("loading muscle\n");
        std::string muscleFileName = object["Muscle filename"].GetString();
        MakeMuscles(std::string(SOFTCON_DIR)+"/data/"+muscleFileName,1.0);
        std::cout << mMuscles.size() << " mMuscles\n";
    }

    if (object.HasMember("Attach Verts")) {
        printf("loading Attach y\n");
        mIsAttach = true;
        mAttachStiffness = object["Attach Stiffness"].GetDouble();
        const auto &verts = object["Attach Verts"];
        for (int i = 0; i < verts.Size(); ++i)
            mAttachIdx.insert(verts[i].GetInt());
    }

    if (object.HasMember("Sampling filename")) {
        printf("loading sampling\n");
        std::string samplingFileName = object["Sampling filename"].GetString();
        std::ifstream ifs(std::string(SOFTCON_DIR)+"/data/"+samplingFileName);
        std::string str;
        std::stringstream ss;
        int index;
        
        str.clear();
        ss.clear();
        std::getline(ifs,str);
        ss.str(str);
        while(ss>>index) {
            mSamplingIndex.push_back(index);
        }
        ifs.close();
    }

    const auto& vertices = mMesh->GetVertices();
    printf("----------- set pneumatic\n");
    mVolumeVertIdx = std::set<int>();
    const auto &volverts = object["Volume Verts"];
    int num_volverts = volverts.Size();
    for (int i = 0; i < num_volverts; ++i) {
        const auto &lines = volverts[i];
        int num_verts = lines.Size();
        std::vector<int> tmp;
        for (int j = 0; j < num_verts; ++j)
            tmp.push_back(lines[j].GetInt());
        for (int j = 0; j < num_verts-1; ++j) {
            int i0 = tmp[j], i1 = tmp[j+1];
            Eigen::Matrix<TinyScalar, 3, 1> p0 = vertices[i0], p1 = vertices[i1];
            TinyScalar length = (p0 - p1).norm();

            auto fem_constraint = new PneumaticConstraint<TinyScalar, TinyConstants>(i0, i1, length, mMuscleStiffness);
            mVolumeMuscles.push_back(fem_constraint);
            mConstraints.push_back(fem_constraint);
        }
        // mVolumeVertIdx.insert(volverts[i].GetInt());
    }


    printf("----------- set rigid\n");
    //set rigids -------------------------------------------
    std::vector<std::vector<int>> rigidBodies;
    {
        const auto &bones = object["Bones"];
        int num_bones = bones.Size();
        for (int i = 0; i < num_bones; ++i) {
            const auto &bone = bones[i];
            int pid = bone["Parent ID"].GetInt();
            std::string jointType = bone["Joint Type"].GetString();
            auto jt = jointType=="ALL" ? Link<TinyScalar, TinyConstants>::ALL : 
                (jointType=="ROTATIONAL" ? Link<TinyScalar, TinyConstants>::ROTATIONAL : 
                (jointType=="PRISMATIC" ? Link<TinyScalar, TinyConstants>::PRISMATIC : 
                    Link<TinyScalar, TinyConstants>::ERR));
            Eigen::Matrix<TinyScalar, 3, 1> jointPos, jointAxis;
            for (int i = 0; i < 3; ++i)
                jointPos[i] = bone["Joint Pos"][i].GetDouble();
            if (pid >= 0)
                jointPos += offset;
            for (int i = 0; i < 3; ++i)
                jointAxis[i] = bone["Joint Axis"][i].GetDouble();
            jointAxis = jointAxis / jointAxis.norm();
            std::vector<int> vertidx;
            const auto &verts = bone["Vertices"];
            for (int i = 0; i < verts.Size(); ++i)
                vertidx.push_back(verts[i].GetInt());

            rigidBodies.push_back(vertidx);
            Link<TinyScalar, TinyConstants> *pa = pid < 0 ? NULL : mlinks[pid];
            mlinks.push_back(new Link<TinyScalar, TinyConstants>(pa, jointPos, jt, jointAxis, i));
            Link<TinyScalar, TinyConstants> *cur = mlinks.back();
            for (int ind : vertidx)
                cur->restpos.push_back(vertices[ind]);
        }
    }
    //set rigids -------------------------------------------


    mRigidBodies = rigidBodies;
    for (auto &idxs : rigidBodies)
        for (int idx : idxs) {
            mrigid.push_back(idx*3+0);
            mrigid.push_back(idx*3+1);
            mrigid.push_back(idx*3+2);
        }
    std::vector<int> mrigid_ori = mrigid;
    std::sort(mrigid_ori.begin(), mrigid_ori.end());
    for (int i = 0, j = 0; i < 3*vertices.size(); ++i) {
        if (j == mrigid_ori.size() || i < mrigid_ori[j])
            mnonrigid.push_back(i);
        if (j < mrigid_ori.size() && i == mrigid_ori[j])
            ++j;
    }

    mNumUnknown = 3*vertices.size();

    mCenterIndex = object["Center Index"].GetInt();
    mEndEffectorIndex = object["End Effector Index"].GetInt();
    for (int i = 0; i < 3; ++i)
        mLocalContourIndex[i] = object["Local Contour Index"][i].GetInt();
}

template <typename TinyScalar, typename TinyConstants> 
void Dfobj<TinyScalar, TinyConstants>::
SetActions(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& actions)
{
    if (mActuatorType == 1) { // skeleton
        mSoftWorld->mJointTorque = actions;
        // std::cout <<"??? actions "<<actions[0]<<" "<<actions[1]<<" "<<actions[2]<<" "<<actions[3]<<std::endl;
    } else if (mActuatorType == 2) { // muscle
        Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> real_actions = mNormalizer->NormToReal(actions);
        SetActivationLevels(real_actions, mPhase);
    } else if (mActuatorType == 3) { // pneumatic
        for (int i = 0; i < mVolumeMuscles.size(); i++) {
            mVolumeMuscles[i]->SetActivationLevel(actions[i]);
        }
    }
}

template <typename TinyScalar, typename TinyConstants> 
void Dfobj<TinyScalar, TinyConstants>::
InitializeActions()
{
    // 0 no actuator
    // 1 skeleton
    // 2 muscle
    // 3 pneumatic
    int num_action;
    if (mActuatorType == 0) {
        printf("no actuator\n");
    } else if (mActuatorType == 1) {
        printf("skeleton actuator\n");
        num_action = mRigidBodies.size() - 1;
        mActions.resize(num_action);
        std::cout << "num_action " << num_action << "\n";
    } else if (mActuatorType == 2) {
        printf("muscle actuator\n");
        const auto& muscles = GetMuscles();
        num_action = 4 * muscles.size();
        mActions.resize(num_action);

        Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> real_lower_bound(num_action);
        Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> real_upper_bound(num_action);

        std::cout << num_action << " num_action\n";
        std::cout << muscles.size() << " muscles\n";
        int cnt = 0;
        for(const auto& m : muscles) 
        {
            real_lower_bound.segment(cnt,4) = m->GetActionLowerBound();
            real_upper_bound.segment(cnt,4) = m->GetActionUpperBound();
            cnt += 4;
        }

        mNormLowerBound.resize(real_lower_bound.size());
        mNormUpperBound.resize(real_upper_bound.size());

        mNormLowerBound.setOnes();
        mNormLowerBound *= -5.0;
        mNormUpperBound.setOnes();
        mNormUpperBound *= 5.0;

        mNormalizer = new Normalizer<TinyScalar, TinyConstants>(real_upper_bound,real_lower_bound,
            mNormUpperBound,mNormLowerBound);
    } else if (mActuatorType == 3) {
        num_action = mVolumeMuscles.size();
        mActions.resize(num_action);
        printf("pneumatic actuator\n");
    }

}

template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, 3, 3> Dfobj<TinyScalar, TinyConstants>::
GetReferenceRotation(bool type,const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x)
{
    int idx0 = mLocalContourIndex[0];
    int idx1 = mLocalContourIndex[1];
    int idx2 = mLocalContourIndex[2];

    Eigen::Matrix<TinyScalar, 3, 1> p0, p1, p2;
    p0 = x.template block<3,1>(3*idx0,0);
    p1 = x.template block<3,1>(3*idx1,0);
    p2 = x.template block<3,1>(3*idx2,0);

    Eigen::Matrix<TinyScalar, 3, 1> v_x, v_y, v_z;
    v_x = (p1-p0).normalized();
    v_y = (p1-p0).cross((p2-p0)).normalized();
    v_z = v_x.cross(v_y).normalized();

    
    Eigen::Matrix<TinyScalar, 3, 3> R;
    R(0,0) = v_x[0];  R(0,1) = v_y[0];  R(0,2) = v_z[0];
    R(1,0) = v_x[1];  R(1,1) = v_y[1];  R(1,2) = v_z[1];
    R(2,0) = v_x[2];  R(2,1) = v_y[2];  R(2,2) = v_z[2];

    if (type == LOCAL)
        return R.transpose();
    else
        return R;
}


template <typename TinyScalar, typename TinyConstants> 
void Dfobj<TinyScalar, TinyConstants>::
SetInitReferenceRotation(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x)
{
    mInitReferenceMatrix = GetReferenceRotation(GLOBAL,x);
}

template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, 3, 1> Dfobj<TinyScalar, TinyConstants>::
GetUpVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x)
{
    Eigen::Matrix<TinyScalar, 3, 3> R = GetReferenceRotation(GLOBAL,x);
    auto v_up_simul = R*mInitReferenceMatrix.inverse().col(1);
    return v_up_simul;
}

template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, 3, 1> Dfobj<TinyScalar, TinyConstants>::
GetForwardVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x)
{
    Eigen::Matrix<TinyScalar, 3, 3> R = GetReferenceRotation(GLOBAL,x);
    return R.col(1);
}


template <typename TinyScalar, typename TinyConstants> 
void Dfobj<TinyScalar, TinyConstants>::
Initialize(FEM::World<TinyScalar, TinyConstants>* world)
{
    const auto& vertices = mMesh->GetVertices();
    const auto& tetrahedras = mMesh->GetTetrahedrons();

    std::vector<Eigen::Vector3i> triangles;
    std::vector<std::pair<Eigen::Vector3i,int>> surfaces;

    m_masses.resize(vertices.size());
    TinyScalar totvol = 0;

    for(const auto& tet : tetrahedras)
    {
        int i0,i1,i2,i3;
        Eigen::Matrix<TinyScalar, 3, 1> p0,p1,p2,p3;
        
        i0 = tet[0];
        i1 = tet[1];
        i2 = tet[2];
        i3 = tet[3];
        p0 = vertices[i0];
        p1 = vertices[i1];
        p2 = vertices[i2];
        p3 = vertices[i3];

        Eigen::Matrix<TinyScalar, 3, 3> Dm;

        Dm.template block<3,1>(0,0) = p1 - p0;
        Dm.template block<3,1>(0,1) = p2 - p0;
        Dm.template block<3,1>(0,2) = p3 - p0;
        if(Dm.determinant()<0)
        {
            i0 = tet[0];
            i1 = tet[2];
            i2 = tet[1];
            i3 = tet[3];
            p0 = vertices[i0];
            p1 = vertices[i1];
            p2 = vertices[i2];
            p3 = vertices[i3];
            Dm.template block<3,1>(0,0) = p1 - p0;
            Dm.template block<3,1>(0,1) = p2 - p0;
            Dm.template block<3,1>(0,2) = p3 - p0;
        }
        bool is_actuate = false;
        // is_actuate = is_actuate || (mVolumeVertIdx.find(i0) != mVolumeVertIdx.end());
        // is_actuate = is_actuate || (mVolumeVertIdx.find(i1) != mVolumeVertIdx.end());
        // is_actuate = is_actuate || (mVolumeVertIdx.find(i2) != mVolumeVertIdx.end());
        // is_actuate = is_actuate || (mVolumeVertIdx.find(i3) != mVolumeVertIdx.end());
        // bool is_actuate =  mIsVolumeAct & p0[1] < -0.1;
        if (is_actuate) {
            // auto fem_constraint = new CorotateFEMConstraintActuate<TinyScalar, TinyConstants>(mYoungsModulus,mPoissonRatio,i0,i1,i2,i3,
            //     1.0/6.0*(Dm.determinant()),Dm.inverse());
            // mVolumeMuscles.push_back(fem_constraint);
            // mConstraints.push_back(fem_constraint);
        }
        else {
            auto fem_constraint = new CorotateFEMConstraint<TinyScalar, TinyConstants>(mYoungsModulus,mPoissonRatio,i0,i1,i2,i3,
                1.0/6.0*(Dm.determinant()),Dm.inverse());
            mConstraints.push_back(fem_constraint);
        }

        if (mIsMassSame) {
            TinyScalar vol = 1.0 / 6.0 * (Dm.determinant());
            totvol += vol;
            for (int i = 0; i < 4; ++i) {
                m_masses[tet[i]] += 0.25 * vol;
            }
        }

        int sorted_idx[4] ={i0,i1,i2,i3};
        std::sort(sorted_idx,sorted_idx+4);
        triangles.push_back(Eigen::Vector3i(sorted_idx[0],sorted_idx[1],sorted_idx[2]));
        triangles.push_back(Eigen::Vector3i(sorted_idx[0],sorted_idx[1],sorted_idx[3]));
        triangles.push_back(Eigen::Vector3i(sorted_idx[0],sorted_idx[2],sorted_idx[3]));
        triangles.push_back(Eigen::Vector3i(sorted_idx[1],sorted_idx[2],sorted_idx[3]));
    }

    if (mIsMassSame) {
        for (std::size_t vId = 0; vId < m_masses.size(); ++vId)
            m_masses[vId] = mTotalMass / totvol * m_masses[vId];
    }
    else {
        for (std::size_t vId = 0; vId < m_masses.size(); ++vId)
            m_masses[vId] = mTotalMass / m_masses.size();
    }


    if (mIsAttach) {
        for (int i = 0; i < vertices.size(); i++) {
            if (mAttachIdx.find(i) != mAttachIdx.end()) {
                auto attach_constraint = new AttachmentConstraint<TinyScalar, TinyConstants>(
                    mAttachStiffness, i, vertices[i]);
                mConstraints.push_back(attach_constraint);
            }
        }
    }

    for(int i=0;i<triangles.size();i++)
    {
        Eigen::Vector3i t_i = triangles[i];
        bool unique = true;
        for(int j=0;j<triangles.size();j++)
        {
            if(i==j)
                continue;
            if(t_i.isApprox(triangles[j]))
                unique = false;
        }
        if(unique)
            surfaces.push_back(std::make_pair(t_i,i/4));
    }
    
    for(int i=0;i<surfaces.size();i++) {
        mContactFace.push_back(surfaces[i].first);
        for (int j=0;j<3;j++) {
            int idx_f = mContactFace[i][j];
            if (std::find(mContactIdx.begin(), mContactIdx.end(), idx_f) 
                    == mContactIdx.end() )
                mContactIdx.push_back(idx_f);
        }
    }
    
    


    for(int i=0;i<surfaces.size();i++)
    {
        Eigen::Vector3i t_i = surfaces[i].first;
        int tet_index = surfaces[i].second;
        
        int i0 = tetrahedras[tet_index][0], i1 = tetrahedras[tet_index][1], i2 = tetrahedras[tet_index][2], i3 = tetrahedras[tet_index][3];
        Eigen::Matrix<TinyScalar, 3, 1> p0 = vertices[i0],p1 = vertices[i1],p2 = vertices[i2],p3 = vertices[i3];
        Eigen::Matrix<TinyScalar, 3, 1> center = 0.25*(p0+p1+p2+p3);

        Eigen::Matrix<TinyScalar, 3, 1> q0 = vertices[t_i[0]],q1 = vertices[t_i[1]],q2 = vertices[t_i[2]];
        Eigen::Matrix<TinyScalar, 3, 1> n = (q1-q0).cross(q2-q1);

        if((center-q0).dot(n)>0.0)
        {
            int j1 = t_i[0];
            t_i[0] = t_i[1];
            t_i[1] = j1;
        }

        mContours.push_back(t_i);
    }
    std::cout << mContours.size() << " contours size \n";
    std::cout << surfaces.size() << " surfaces size \n";

    for(const auto& m : mMuscles)
        for(const auto& s : m->GetSegments()) 
            for(int i=0; i<s->GetMuscleConstraints().size(); i++)
                mConstraints.push_back(s->GetMuscleConstraints()[i]);

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> v(vertices.size()*3);
    for(int i =0;i<vertices.size();i++)
        v.template block<3,1>(i*3,0) = vertices[i];


    world->AddBody(v,mConstraints,m_masses,mContactIdx,mContactFace, 
    	mRigidBodies, mnonrigid, mrigid, mlinks);
    if (mMuscles.size() > 0) {
        SetKey("Dfobj.param");
    }

    mSoftWorld = world;
}

template <typename TinyScalar, typename TinyConstants> 
void Dfobj<TinyScalar, TinyConstants>::
OutputSurface(std::string filename, const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &x) {
    FILE *file = fopen(filename.c_str(), "w");
    const auto& vertices = mMesh->GetVertices();
    for (int i = 0; i < vertices.size(); ++i)
        fprintf(file, "v %06f %06f %06f\n",
            TinyConstants::getDouble(x[i*3+0]),
            TinyConstants::getDouble(x[i*3+1]),
            TinyConstants::getDouble(x[i*3+2]));
    for (int i = 0; i < mContours.size(); ++i)
        fprintf(file, "f %d %d %d\n",mContours[i][0]+1,mContours[i][1]+1,mContours[i][2]+1);
    fclose(file);
}

template <typename TinyScalar, typename TinyConstants> 
void Dfobj<TinyScalar, TinyConstants>::
SetKey(std::string filename)
{
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> mGivenKey,key;

    int row=50, col=mMuscles.size();
    
    mGivenKey.resize(row,col);

    for(int i=0; i<20; i++) {
        mGivenKey.row(i)<<0.025,-0.0625,0.025,-0.0625,0.025,-0.0625,0.025,-0.0625,0.025,-0.0625,0.025,-0.0625,0.025,-0.0625,0.025,-0.0625;
    }
    for(int i=20; i<30; i++) {
        mGivenKey.row(i)<<-0.125,0.03125,-0.125,0.03125,-0.125,0.03125,-0.125,0.03125,-0.125,0.03125,-0.125,0.03125,-0.125,0.03125,-0.125,0.03125;
    }
    for(int i=30; i<50; i++) {
        mGivenKey.row(i)<<0.0,-0.125,0.0,-0.125,0.0,-0.125,0.0,-0.125,0.0,-0.125,0.0,-0.125,0.0,-0.125,0.0,-0.125;
    }

    key.resize(mGivenKey.rows(),mGivenKey.cols());
    key.setZero();

    for(int j=0; j<key.cols(); j++) {
        for(int i=0; i<key.rows(); i++) {
            if(i==0) {
                key(i,j) += mGivenKey(i,j);
            } else {
                key(i,j) += key(i-1,j) + mGivenKey(i,j);
            }

            key(i,j) = std::min(std::max(key(i,j),TinyConstants::zero()),TinyConstants::one());
        }
    }

    for(int i=0; i<key.cols(); i++) {
        mMuscles[i]->SetKey(key.col(i));
    }  
    
}


template <typename TinyScalar, typename TinyConstants> 
void Dfobj<TinyScalar, TinyConstants>::
SetActivationLevels(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& actions,const int& phase)
{
    int idx = 0;
    // std::cout << actions[0] << " " << actions[1] << " " << actions[2] << " " << actions[3] << "\n";
    for(int i = 0; i < mMuscles.size(); i++)
    {
        mMuscles[i]->SetActivationLevels(actions.segment(idx,4), phase);
        idx += 4;
    }
    // static int  count = 0;
    // TinyScalar action = 1 + std::max(1000*sin(count++ / 50.), 0.);
    // if (mIsVolumeAct)
    //      for (int i=0; i<mVolumeMuscles.size(); i++)
    //         mVolumeMuscles[i]->SetActivationLevel(action);
}

template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> 
Dfobj<TinyScalar, TinyConstants>::
GetActivationLevelsAggregate()
{
    int len = 0;
    int idx = 0;
    for(int i=0; i<mMuscles.size(); i++)
        len += mMuscles[i]->mActivationLevels.size();

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ans(len);
    for(int i=0; i<mMuscles.size(); i++) {
        const auto& a = mMuscles[i]->mActivationLevels;
        for (int j = 0; j < a.size(); j++)
            ans[idx + j] = a[j];
        idx += a.size();
    }
    return ans;
}

template <typename TinyScalar, typename TinyConstants> 
void Dfobj<TinyScalar, TinyConstants>::
SetActivationLevelsAggregate(Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& v)
{
    int idx = 0;
    for(int i=0; i<mMuscles.size(); i++) {
        auto& a = mMuscles[i]->mActivationLevels;
        for (int j = 0; j < a.size(); j++)
            a[j] = v[idx + j];
        idx += a.size();
    }
}


template <typename TinyScalar, typename TinyConstants> 
void Dfobj<TinyScalar, TinyConstants>::
SetMesh(const std::string& path,const Eigen::Transform<TinyScalar, 3, 2>& T)
{
    std::cout << std::string(SOFTCON_DIR)+"/data/"+path,mScalingMatrix;
    mScalingMatrix = T;
    mMesh = new OBJMesh<TinyScalar, TinyConstants>(std::string(SOFTCON_DIR)+"/data/"+path,mScalingMatrix);
}


template <typename TinyScalar, typename TinyConstants> 
void Dfobj<TinyScalar, TinyConstants>::
MakeMuscles(const std::string& path,const TinyScalar& gamma)
{
    printf("MakeMuscle! \n");
    std::ifstream ifs(path);
    if(!ifs.is_open()){
        std::cout << "Muscle file doesn't exist." << std::endl;
        exit(0);
    }
    std::string str;
    std::string index;
    std::stringstream ss;
    
    int fiber_type;
    int num_sampling;

    int current_fiber_index = 0;
    
    std::vector<Eigen::Matrix<TinyScalar, 3, 1>> point_list;

    while(!ifs.eof())
    {
        str.clear();
        index.clear();
        ss.clear();

        std::getline(ifs,str);
        ss.str(str);

        int fiber_index;
        ss>>fiber_index;

        if(fiber_index != current_fiber_index) {
            fiber_type = 0;
            num_sampling = 50;

            mMuscles.push_back(new Muscle<TinyScalar, TinyConstants>(num_sampling,point_list));
            Muscle<TinyScalar, TinyConstants>* muscle = mMuscles.back();
            muscle->Initialize(mMesh,mMuscleStiffness);
            point_list.clear();
            current_fiber_index = fiber_index;
        }

        int segment_index;
        ss>>segment_index;

        Eigen::Matrix<TinyScalar, 3, 1> point;
        ss>>point[0]>>point[1]>>point[2];

        point = mScalingMatrix*point;
        point_list.push_back(point);
    }
    ifs.close();

    mMuscles.push_back(new Muscle<TinyScalar, TinyConstants>(num_sampling,point_list));
    Muscle<TinyScalar, TinyConstants>* muscle = mMuscles.back();
    muscle->Initialize(mMesh,mMuscleStiffness);
}
template <typename TinyScalar, typename TinyConstants> 
TinyScalar Dfobj<TinyScalar, TinyConstants>::
C_d(TinyScalar theta)
{
    return -cos(theta*2.0)+1.05;
}
template <typename TinyScalar, typename TinyConstants> 
TinyScalar Dfobj<TinyScalar, TinyConstants>::
C_t(TinyScalar theta)
{
    return 0.25*(exp(0.8*theta) -1);
}
template <typename TinyScalar, typename TinyConstants> 
TinyScalar Dfobj<TinyScalar, TinyConstants>::
Cd(TinyScalar aoa,TinyScalar a,TinyScalar c) {
    return (a-c) * (0.5*(-cos(2.0*aoa)+1.0)) + c;
}
template <typename TinyScalar, typename TinyConstants> 
TinyScalar Dfobj<TinyScalar, TinyConstants>::
_interpolate(TinyScalar t, TinyScalar a, TinyScalar b)
{
    return  (1.0 - t) * a + t * b;
}
template <typename TinyScalar, typename TinyConstants> 
TinyScalar Dfobj<TinyScalar, TinyConstants>::
Cl(TinyScalar aoa, TinyScalar cutoff, TinyScalar x[5], TinyScalar y[5]) {
    const   TinyScalar  xa = x[0],   ya = y[0];
    const   TinyScalar  xb = x[1],   yb = y[1];
    const   TinyScalar  xc = x[2],   yc = y[2];
    const   TinyScalar  xd = x[3],   yd = y[3];
    const   TinyScalar  xe = x[4],   ye = y[4];

    TinyScalar  theta = aoa * 180.0 / M_PI;

    if (fabs(theta) > cutoff)
        return  0.0;

    if (xa <= theta && theta < xb)
        return _interpolate((theta - xa) / (xb - xa), ya, yb);
    else if (xb <= theta && theta < xc)
        return _interpolate((theta - xb) / (xc - xb), yb, yc);
    else if (xc <= theta && theta <= xd)
        return _interpolate((theta - xc) / (xd - xc), yc, yd);
    else if (xd <= theta && theta <= xe)
        return _interpolate((theta - xd) / (xe - xd), yd, ye);
    else
    {
        std::cout << "Error: this should not be reached... " << std::endl;
        std::cout << "Theta: " << aoa << "(deg) " << theta << "(rad)" << std::endl;
        std::cout << "x: " << xa << " " << xb << " " << xc << " " << xd << std::endl;
        std::cout << "y: " << ya << " " << yb << " " << yc << " " << yd << std::endl;
        std::cin.get();
    }
    return 0.0;
}

template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> Dfobj<TinyScalar, TinyConstants>::
ComputeDragForces(FEM::World<TinyScalar, TinyConstants>* world)
{
    // const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x = world->GetPositions();
    // const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x_dot = world->GetVelocities();

    // int n = x.rows();
    // Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> f = Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>::Zero(n);
    // Eigen::Matrix<TinyScalar, 3, 1> f_sum = Eigen::Matrix<TinyScalar, 3, 1>::Zero(); 
    // Eigen::Matrix<TinyScalar, 3, 1> avg_vel = Eigen::Matrix<TinyScalar, 3, 1>::Zero();

    // for(const auto& con : mContours)
    // {
    //  int i0 = con[0];
    //  int i1 = con[1];
    //  int i2 = con[2];

    //  Eigen::Matrix<TinyScalar, 3, 1> p0 = x.segment<3>(i0*3);
    //  Eigen::Matrix<TinyScalar, 3, 1> p1 = x.segment<3>(i1*3);
    //  Eigen::Matrix<TinyScalar, 3, 1> p2 = x.segment<3>(i2*3);
    //  Eigen::Matrix<TinyScalar, 3, 1> com = (p0+p1+p2)/3.0;

    //  Eigen::Matrix<TinyScalar, 3, 1> v0 = x_dot.segment<3>(i0*3);
    //  Eigen::Matrix<TinyScalar, 3, 1> v1 = x_dot.segment<3>(i1*3);
    //  Eigen::Matrix<TinyScalar, 3, 1> v2 = x_dot.segment<3>(i2*3);

    //  Eigen::Matrix<TinyScalar, 3, 1> v = -(v0+v1+v2)/3.0;
    //  if(v.norm()<1E-6)
    //      continue;
    //  Eigen::Matrix<TinyScalar, 3, 1> n = (p1-p0).cross(p2-p1);
        
    //  TinyScalar area = 0.5*n.norm();
    //  n.normalize();
    //  n = -n;
    
    //  TinyScalar theta = atan2(v.dot(n),(v-v.dot(n)*n).norm());
    //  Eigen::Matrix<TinyScalar, 3, 1> d = v.normalized();
    //  Eigen::Matrix<TinyScalar, 3, 1> t = -n;
    //  Eigen::Matrix<TinyScalar, 3, 1> fv = GetForwardVector(x);
    //  TinyScalar theta2 = atan2(t.dot(fv),(t-t.dot(fv)*fv).norm());
    //  Eigen::Matrix<TinyScalar, 3, 1> l = d.cross(n);

    //  l.normalize();
    //  l = l.cross(d);

    //  TinyScalar f_d = 0.5*1000.0*area*v.squaredNorm()*C_d(theta)*0.5;
    //  TinyScalar f_t = 0.5*1000.0*area*v.squaredNorm()*C_t(theta2)*std::abs(fv.dot(t))*2;
    //  Eigen::Matrix<TinyScalar, 3, 1> f_i = 0.333*(f_d*d+f_t*t);
    //  f.segment<3>(i0*3) += f_i;
    //  f.segment<3>(i1*3) += f_i;
    //  f.segment<3>(i2*3) += f_i;
    // }
    // TinyScalar clip_f = 1E3;
    // for(int i =0;i<n;i++)
    // {
    //  f[i] = std::max(-clip_f,std::min(clip_f,f[i]));
    // }
    // return f;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> force_drag = Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>::Zero(3*world->GetNumVertices());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> force_lift = Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>::Zero(3*world->GetNumVertices());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> force_total = Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>::Zero(3*world->GetNumVertices());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> velocities = world->GetVelocities();
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> positions = world->GetPositions();

    // if (mIsUnderwater) {
        TinyScalar max_force = 1.0e03;

        TinyScalar cl_x[5] = {-90, -10, -5, 15, 90};
        TinyScalar cl_y[5] = {-0.5, -1.26, 0.0, 1.8, 0.5};
        TinyScalar cutoff = 85.0;

        for(int i=0; i<mContours.size(); i++)
        {
            Eigen::Matrix<TinyScalar, 3, 1> p1 = positions.template segment<3>(mContours[i][0]*3);
            Eigen::Matrix<TinyScalar, 3, 1> p2 = positions.template segment<3>(mContours[i][1]*3);
            Eigen::Matrix<TinyScalar, 3, 1> p3 = positions.template segment<3>(mContours[i][2]*3);
            Eigen::Matrix<TinyScalar, 3, 1> v1 = velocities.template segment<3>(mContours[i][0]*3);
            Eigen::Matrix<TinyScalar, 3, 1> v2 = velocities.template segment<3>(mContours[i][1]*3);
            Eigen::Matrix<TinyScalar, 3, 1> v3 = velocities.template segment<3>(mContours[i][2]*3);

            Eigen::Matrix<TinyScalar, 3, 1> v_water(0.0, 0.0, 0.0);
            Eigen::Matrix<TinyScalar, 3, 1> v_rel = v_water - (v1+v2+v3)/3.0;
            if(v_rel.norm() < 1e-6) continue;
            Eigen::Matrix<TinyScalar, 3, 1> v_rel_norm = v_rel.normalized();

            Eigen::Matrix<TinyScalar, 3, 1> n = ((p2-p1).cross(p3-p1)).normalized();
            Eigen::Matrix<TinyScalar, 3, 1> d = -n;
            Eigen::Matrix<TinyScalar, 3, 1> d_lift(0.0, 0.0, 0.0);

            Eigen::Matrix<TinyScalar, 3, 1> force_drag_per_face(0.0, 0.0, 0.0);
            Eigen::Matrix<TinyScalar, 3, 1> force_lift_per_face(0.0, 0.0, 0.0);

            TinyScalar area = 0.5*((p2-p1).cross(p3-p1)).norm();

            if (area < 1.0e-10)
            {
                std::cerr << "Error: Area is too small, you should check the simulation" << std::endl;
                std::cin.get();
            }

            TinyScalar d_dot_v_rel_norm = d.dot(v_rel_norm);
            TinyScalar aoa = 0.5 * M_PI - acos(d_dot_v_rel_norm);
            TinyScalar f = 1.0;
        
            // Ignore faces which have reverse direction 
            if (d_dot_v_rel_norm > 0.0)
            {
                if (d_dot_v_rel_norm > 1.0 || d_dot_v_rel_norm < -1.0)
                {
                    std::cerr << "Error: this should not be reached... " << std::endl;
                    // std::cerr << "d_dot_v_rel_norm: " << d_dot_v_rel_norm << std::endl;
                    // std::cerr << "Drag: " << d.transpose() << std::endl;
                    // std::cerr << "Vel: " << v_rel_norm.transpose() << std::endl;
                    std::cin.get();
                } 

                if (fabs(aoa) > 0.5 * M_PI)
                {
                    std::cerr << "Error: this should not be reached... " << std::endl;
                    // std::cerr << "aoa: " << aoa << std::endl;
                    std::cin.get();
                }

                force_drag_per_face = 0.5 * 1000.0 * (1.0*area) * Cd(aoa,f,0.0) * v_rel.squaredNorm() * v_rel_norm;
            }

            if (d_dot_v_rel_norm > 0.0)
            {
                d_lift = (v_rel_norm.cross(d)).normalized().cross(v_rel_norm);
                force_lift_per_face = 0.5 * 1000.0 * (1.0*area) * Cl(aoa,cutoff,cl_x,cl_y) * v_rel.squaredNorm() * d_lift;
            }

            if (force_drag_per_face.norm() > max_force)
            {
                // std::cerr << "Error: max_force reached..." << std::endl;
                force_drag_per_face = max_force * force_drag_per_face.normalized();
            }

            if (force_lift_per_face.norm() > max_force)
            {
                // std::cerr << "Error: max_force reached..." << std::endl;
                force_lift_per_face = max_force * force_lift_per_face.normalized();
            }

            force_total.template segment<3>(mContours[i][0]*3) += force_drag_per_face/3.0 + force_lift_per_face/3.0;
            force_total.template segment<3>(mContours[i][1]*3) += force_drag_per_face/3.0 + force_lift_per_face/3.0;
            force_total.template segment<3>(mContours[i][2]*3) += force_drag_per_face/3.0 + force_lift_per_face/3.0;
        }
    // }
    return force_total;
}


template <typename TinyScalar, typename TinyConstants> 
void Dfobj<TinyScalar, TinyConstants>::
SetVertexNormal()
{
    const auto& vertices = mMesh->GetVertices();
    const auto& tetrahedras = mMesh->GetTetrahedrons();

    mVertexNormal.resize(3*vertices.size());
    mVertexNormal.setZero();
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> cnt_list(vertices.size());
    cnt_list.setZero();

    for(int i=0;i<mContours.size();i++)
    {
        int i0 = mContours[i][0];
        int i1 = mContours[i][1];
        int i2 = mContours[i][2];

        Eigen::Matrix<TinyScalar, 3, 1> p0 = vertices[i0];
        Eigen::Matrix<TinyScalar, 3, 1> p1 = vertices[i1];
        Eigen::Matrix<TinyScalar, 3, 1> p2 = vertices[i2];

        Eigen::Matrix<TinyScalar, 3, 1> face_normal;
        face_normal = (p1-p0).cross(p2-p1);

        mVertexNormal.template block<3,1>(3*i0,0) = cnt_list[i0]/(cnt_list[i0]+1)*mVertexNormal.template block<3,1>(3*i0,0) + 1/(cnt_list[i0]+1)*face_normal;
        mVertexNormal.template block<3,1>(3*i1,0) = cnt_list[i1]/(cnt_list[i1]+1)*mVertexNormal.template block<3,1>(3*i1,0) + 1/(cnt_list[i1]+1)*face_normal;
        mVertexNormal.template block<3,1>(3*i2,0) = cnt_list[i2]/(cnt_list[i2]+1)*mVertexNormal.template block<3,1>(3*i2,0) + 1/(cnt_list[i2]+1)*face_normal;

        cnt_list[i0] += 1;
        cnt_list[i1] += 1;
        cnt_list[i2] += 1;
    }
}
#endif