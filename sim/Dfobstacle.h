#ifndef __Dfobstacle_H__
#define __Dfobstacle_H__
#include "Dfobj.h"

#include <fstream>
#include <math.h>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>


#include "utils/jsonHelper.h"
#include "utils/geometry_utils.h"

template <typename TinyScalar, typename TinyConstants> 
class Dfobstacle : public Dfobj<TinyScalar, TinyConstants>
{
public:
    using SurfaceMesh = CGAL::Surface_mesh<Eigen::Matrix<TinyScalar, 3, 1>>;

    Dfobstacle();
    ~Dfobstacle();

    void Initialize(FEM::World<TinyScalar, TinyConstants>* world);
    void SetDfobj(const rapidjson::Value &object, Eigen::Matrix<TinyScalar, 3, 1> offset);
    void SetMesh(const std::string& path,const  Eigen::Transform<TinyScalar, 3, 2>& T);
    void MakeMuscles(const std::string& path,const double& gamma);

    void SetActivationLevels(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& actions,const int& phase); 
    void SetKey(std::string filename);

    Eigen::Matrix<TinyScalar, 3, 3> GetReferenceRotation(bool type,const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
    void SetInitReferenceRotation(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);

    const std::vector<Eigen::Vector3i>& GetContours(){return mContours;};

    void SetVertexNormal();
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& GetVertexNormal(){return mVertexNormal;};

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ComputeDragForces(FEM::World<TinyScalar, TinyConstants>* world);

    Eigen::Matrix<TinyScalar, 3, 1> GetUpVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
    Eigen::Matrix<TinyScalar, 3, 1> GetForwardVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x); 

    const int& GetCenterIndex() {return mCenterIndex;};
    const int& GetEndEffectorIndex() {return mEndEffectorIndex;};
    const std::vector<int>& GetSamplingIndex() {return mSamplingIndex;};

    std::vector<std::vector<TinyScalar>> getFrictionCoefficientsFromArray(const rapidjson::Value& array);


// private:
    double                                  mMuscleStiffness;
    double                                  mYoungsModulus;
    double                                  mPoissonRatio;

    FEM::Mesh<TinyScalar, TinyConstants>*                               mMesh;
     Eigen::Transform<TinyScalar, 3, 2>                         mScalingMatrix;

    std::vector<FEM::Constraint<TinyScalar, TinyConstants>*>            mConstraints;

    std::vector<Eigen::Vector3i>            mContours;

    int                                     mCenterIndex;
    int                                     mEndEffectorIndex;
    std::vector<int>                        mSamplingIndex;
    Eigen::Matrix<TinyScalar, 3, 1>                         mLocalContourIndex;

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                         mVertexNormal;

    Eigen::Matrix<TinyScalar, 3, 3>                         mInitReferenceMatrix;

    bool                                    mIsVolumeAct;
    std::vector<FEM::CorotateFEMConstraintActuate<TinyScalar, TinyConstants>*> mVolumeMuscles;


};


using namespace FEM;

template <typename TinyScalar, typename TinyConstants> 
Dfobstacle<TinyScalar, TinyConstants>::
Dfobstacle()
    :Dfobj<TinyScalar, TinyConstants>(0,0,0)
{
    printf("construct Dfobstacle\n");
}

template <typename TinyScalar, typename TinyConstants> 
Dfobstacle<TinyScalar, TinyConstants>::
~Dfobstacle() {
    delete mMesh;
}

template <typename TinyScalar, typename TinyConstants> 
void Dfobstacle<TinyScalar, TinyConstants>::
SetDfobj(const rapidjson::Value &object, Eigen::Matrix<TinyScalar, 3, 1> offset)
{
    // std::cout <<"\n SetDfobj Dfobstacle + Meta Data Path: " << path << std::endl;
    // std::FILE* file = std::fopen(path.c_str(), "r");
    // if (file == nullptr)
    //     throw std::invalid_argument(std::string("Couldn't open file ") + path);
   
    // constexpr std::size_t buffer_size = 1024;
    // char buffer[buffer_size];

    // rapidjson::FileReadStream&& is = rapidjson::FileReadStream(file, buffer, buffer_size);
    // rapidjson::Document scene_description;
    // scene_description.ParseStream(is);
    // if (!scene_description.IsObject())
    //     throw std::invalid_argument("JSON reader - Error : root must be an object");


    // const rapidjson::Value& mesh_array = jsonRequireArrayCheck(scene_description, "Meshes");
    
    // const auto& object = mesh_array[0];
    std::string rest_state_file_name = jsonRequireString(object, "Obj filename");


    Eigen::Transform<TinyScalar, 3, 2> T = Eigen::Transform<TinyScalar, 3, 2>::Identity();
    T(0,0) *= 1; T(1,1) *= 1; T(2,2) *= 1; 
    T.translation() = offset;
    SetMesh(rest_state_file_name,T);
    this->mNumUnknown = 0;
}

template <typename TinyScalar, typename TinyConstants> 
void Dfobstacle<TinyScalar, TinyConstants>::
Initialize(FEM::World<TinyScalar, TinyConstants>* world)
{

    const auto& vertices = mMesh->GetVertices();



    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> v(vertices.size()*3);
    for(int i =0;i<vertices.size();i++) {
        v.template block<3,1>(i*3,0) = vertices[i];

    }


    for(int i=0;i<mMesh->mTriangles.size();i++) {
        this->mContactFace.push_back(mMesh->mTriangles[i]);
        for (int j=0;j<3;j++) {
            int idx_f = this->mContactFace[i][j];
            if (std::find(this->mContactIdx.begin(), this->mContactIdx.end(), idx_f) 
                    == this->mContactIdx.end() )
                this->mContactIdx.push_back(idx_f);
        }
    }
    sort(this->mContactIdx.begin(), this->mContactIdx.end());
    assert(this->mContactIdx[0] == 0 && this->mContactIdx.size() == vertices.size());


    world->AddBody(v,this->mContactIdx,this->mContactFace);
}


template <typename TinyScalar, typename TinyConstants> 
void Dfobstacle<TinyScalar, TinyConstants>::
SetMesh(const std::string& path,const  Eigen::Transform<TinyScalar, 3, 2>& T)
{
    mScalingMatrix = T;
    mMesh = new OBJClothMesh<TinyScalar, TinyConstants>(
        std::string(SOFTCON_DIR)+"/data/"+path,mScalingMatrix
    );
    auto& vertices = mMesh->GetVertices();
    for(auto& v : mMesh->mVertices)
        v = T*v;
}


template <typename TinyScalar, typename TinyConstants> 
void Dfobstacle<TinyScalar, TinyConstants>::
SetVertexNormal()
{
    const auto& vertices = mMesh->GetVertices();
    const auto& mTriangles = mMesh->mTriangles;

    mVertexNormal.resize(3*vertices.size());
    mVertexNormal.setZero();
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> cnt_list(vertices.size());
    cnt_list.setZero();

    for(int i=0;i<mTriangles.size();i++)
    {
        int i0 = mTriangles[i][0];
        int i1 = mTriangles[i][1];
        int i2 = mTriangles[i][2];

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