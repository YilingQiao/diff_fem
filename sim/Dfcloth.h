#ifndef __Dfcloth_H__
#define __Dfcloth_H__
#include "Dfobj.h"

#include <fstream>
#include <math.h>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>


#include "utils/jsonHelper.h"
#include "utils/geometry_utils.h"
// #include "Dfobstacle.h"
// #include "fem/ClothMesh/ClothPhysicsMesh.h"

// #define LOCAL    0
// #define GLOBAL   1

template <typename TinyScalar, typename TinyConstants> 
class Dfcloth : public Dfobj<TinyScalar, TinyConstants>
{
public:
    using SurfaceMesh = CGAL::Surface_mesh<Eigen::Matrix<TinyScalar, 3, 1>>;

    Dfcloth(const TinyScalar& muscle_stiffness,const TinyScalar& youngs_modulus,const TinyScalar& poisson_ratio);
    ~Dfcloth();

    void Initialize(FEM::World<TinyScalar, TinyConstants>* world);
    void SetDfobj(const rapidjson::Value &object);
    void SetMesh(const std::string& path,const Eigen::Transform<TinyScalar, 3, 2>& T);
    void MakeMuscles(const std::string& path,const TinyScalar& gamma);

    void SetActivationLevels(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& actions,const int& phase); 
    void SetKey(std::string filename);

    Eigen::Matrix<TinyScalar, 3, 3> GetReferenceRotation(bool type,const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
    void SetInitReferenceRotation(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);

    // const std::vector<Eigen::Vector3i>& GetContours(){return mContours;};

    // void SetVertexNormal();
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& GetVertexNormal(){return mVertexNormal;};

    const std::vector<Muscle<TinyScalar, TinyConstants>*>& GetMuscles() {return mMuscles;};

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ComputeDragForces(FEM::World<TinyScalar, TinyConstants>* world);

    Eigen::Matrix<TinyScalar, 3, 1> GetUpVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);
    Eigen::Matrix<TinyScalar, 3, 1> GetForwardVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);    

    const int& GetCenterIndex() {return mCenterIndex;};
    const int& GetEndEffectorIndex() {return mEndEffectorIndex;};
    const std::vector<int>& GetSamplingIndex() {return mSamplingIndex;};

    std::vector<std::vector<TinyScalar>> getFrictionCoefficientsFromArray(const rapidjson::Value& array);
    // std::vector<ClothPhysicMesh<TinyScalar, TinyConstants>> getPhysicMeshesFromArray(const rapidjson::Value& array);
    // ClothPhysicMesh<TinyScalar, TinyConstants> getPhysicMeshFromObject(const rapidjson::Value& array);


    ClothPhysicParameters<TinyScalar, TinyConstants> mParameter;
// private:
    TinyScalar                                  mMuscleStiffness;
    TinyScalar                                  mYoungsModulus;
    TinyScalar                                  mPoissonRatio;

    // FEM::Mesh<TinyScalar, TinyConstants>*                                mMesh;
    Eigen::Transform<TinyScalar, 3, 2>                          mScalingMatrix;

    std::vector<FEM::Constraint<TinyScalar, TinyConstants>*>            mConstraints;

    std::vector<Muscle<TinyScalar, TinyConstants>*>                 mMuscles;

    // std::vector<Eigen::Vector3i>             mContours;

    int                                     mCenterIndex;
    int                                     mEndEffectorIndex;
    std::vector<int>                        mSamplingIndex;
    Eigen::Matrix<TinyScalar, 3, 1>                         mLocalContourIndex;

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                            mVertexNormal;

    Eigen::Matrix<TinyScalar, 3, 3>                         mInitReferenceMatrix;

    bool                                    mIsVolumeAct;
    std::vector<FEM::CorotateFEMConstraintActuate<TinyScalar, TinyConstants>*> mVolumeMuscles;

    std::vector<int>                        mCGAL_vertices_indices;
    SurfaceMesh                             mCGALmesh;
    CollisionMesh<TinyScalar, TinyConstants>*   mCollisionMesh;
    // using ObstaclePtr = std::unique_ptr<Obstacle>;
    // std::vector<ObstaclePtr>                mObstacles;      

    TinyScalar C_d(TinyScalar theta);
    TinyScalar C_t(TinyScalar theta);
    TinyScalar Cd(TinyScalar aoa,TinyScalar a,TinyScalar c); 
    TinyScalar _interpolate(TinyScalar t, TinyScalar a, TinyScalar b);
    TinyScalar Cl(TinyScalar aoa, TinyScalar cutoff, TinyScalar x[5], TinyScalar y[5]); 


};


using namespace FEM;

template <typename TinyScalar, typename TinyConstants> 
Dfcloth<TinyScalar, TinyConstants>::
Dfcloth(const TinyScalar& muscle_stiffness,const TinyScalar& youngs_modulus,const TinyScalar& poisson_ratio)
    :Dfobj<TinyScalar, TinyConstants>(muscle_stiffness, youngs_modulus, poisson_ratio)
{
    printf("construct Dfcloth\n");
}

template <typename TinyScalar, typename TinyConstants> 
Dfcloth<TinyScalar, TinyConstants>::
~Dfcloth() {
    delete this->mMesh;
    delete mCollisionMesh;
}

template <typename TinyScalar, typename TinyConstants> 
void Dfcloth<TinyScalar, TinyConstants>::
SetDfobj(const rapidjson::Value &object)
{
    // std::cout <<"\n SetDfobj Dfcloth + Meta Data Path: " << path << std::endl;
    // std::FILE* file = std::fopen(path.c_str(), "r");
 //    if (file == nullptr)
 //        throw std::invalid_argument(std::string("Couldn't open file ") + path);
   
 //    constexpr std::size_t buffer_size = 1024;
 //    char buffer[buffer_size];

 //    rapidjson::FileReadStream&& is = rapidjson::FileReadStream(file, buffer, buffer_size);
 //    rapidjson::Document scene_description;
 //    scene_description.ParseStream(is);
 //    if (!scene_description.IsObject())
 //        throw std::invalid_argument("JSON reader - Error : root must be an object");


    const rapidjson::Value& friction_coefficients_v =
      jsonRequireArrayCheck(object, "Friction coefficients");
    auto friction_table = getFrictionCoefficientsFromArray(friction_coefficients_v);
    for (auto i:friction_table)
        for (auto j:i)
            std::cout << j << " friction_coefficients\n";


 //    const rapidjson::Value& mesh_array = jsonRequireArrayCheck(scene_description, "Meshes");
    
 //    const auto& object = mesh_array[0];
    std::string rest_state_file_name = jsonRequireString(object, "Obj filename");
    mParameter.area_density = TinyConstants::scalar_from_double(
        jsonRequireDouble(object, "Area Density"));
    mParameter.bend = TinyConstants::scalar_from_double(
        jsonRequireDouble(object, "Bending"));
    mParameter.stretch = TinyConstants::scalar_from_double(
        jsonRequireDouble(object, "Stretching"));
    std::size_t material_identifier = jsonRequireUint(object, "Material Identifier");


    if (object.HasMember("Attach Verts")) {
        printf("loading Attach Verts\n");
        this->mIsAttach = true;
        this->mAttachStiffness = object["Attach Stiffness"].GetDouble();
        const auto &verts = object["Attach Verts"];
        for (int i = 0; i < verts.Size(); ++i)
            this->mAttachIdx.insert(verts[i].GetInt());
    }



    double scale = 0.05;
    if (object.HasMember("Scale"))
        scale = object["Scale"].GetDouble();
    std::cout << scale << " scale\n" ;
    assert(scale > 0);
    Eigen::Transform<TinyScalar, 3, 2> T=Eigen::Transform<TinyScalar, 3, 2>::Identity();
    T(0,0) *= scale; T(1,1) *= scale; T(2,2) *= scale; 
    SetMesh(rest_state_file_name,T);

    const auto& vertices = this->mMesh->GetVertices();
    this->mNumUnknown = 3*vertices.size();
    // for (auto& mesh : getPhysicMeshesFromArray(mesh_array))
    // {
    //     scene.addMesh(std::make_unique<PhysicMesh>(std::move(mesh)));
    // }

    // const rapidjson::Value& obstacle_array = jsonRequireArrayCheck(scene_description, "Obstacles");
    // for (auto& obstacle : getObstaclesFromArray(obstacle_array))
    // {
    //     mObstacles.push_back(std::move(obstacle));
    // }

    // scene.finalize();

    // std::fclose(file);
}

template <typename TinyScalar, typename TinyConstants> 
void Dfcloth<TinyScalar, TinyConstants>::
Initialize(FEM::World<TinyScalar, TinyConstants>* world)
{

    std::cout << mCollisionMesh->m_mesh.number_of_edges() << " m_mesh.number_of_edges() \n";
    std::cout << mCollisionMesh->m_mesh.number_of_vertices() << " m_mesh.number_of_vertices () \n";
    std::cout << mCollisionMesh->m_mesh.number_of_faces() << " m_mesh.number_of_faces () \n";

    std::cout << mConstraints.size() << " mConstraints.size() 1\n";
    mCollisionMesh->initializeClothConstraints(mConstraints);
    std::cout << mConstraints.size() << " mConstraints.size() 2\n";
    // std::cout << edges.size() << " edges.size() \n";
    // for (const auto& e : edges) {
    //     TinyScalar rest_length = mMesh->getInitialEdgeLength(e);
    //     auto edge_constraint = new SpringConstraint<TinyScalar, TinyConstants>(
    //         e[0], e[1], rest_length,  mParameter.stretch);
    //     mConstraints.push_back(edge_constraint);
    // }


    const auto& vertices = this->mMesh->GetVertices();
    // const auto& edges = mMesh->mEdges;

    // // TODO initialize mass
    // // m_masses.resize(vertices.size());
    // // for (std::size_t vId = 0; vId < m_masses.size(); ++vId)
    // // {
    // //     auto triangle_around_vertex_mass_range =
    // //       getTriangleAroundVertexMassRange(vId);
    // //     m_masses[vId] = std::accumulate(boost::begin(triangle_around_vertex_mass_range),
    // //                                     boost::end(triangle_around_vertex_mass_range),
    // //                                     0.) / 3;
    // // }


    // // printf("edges\n");
    // for (const auto& e : edges) {
    //     TinyScalar rest_length = mMesh->getInitialEdgeLength(e);
    //     auto edge_constraint = new SpringConstraint<TinyScalar, TinyConstants>(
    //         e[0], e[1], rest_length,  mParameter.stretch);
    //     mConstraints.push_back(edge_constraint);
    // }

    // for (int i = 0; i < vertices.size(); ++i)
    // {
    //     auto vertices_around_vertex = mMesh->mAdjacentVerts[i];
    //     std::vector<int> vertex_indices;
    //     vertex_indices.push_back(i);
    //     vertex_indices.insert(vertex_indices.end(),
    //                         vertices_around_vertex.begin(),
    //                         vertices_around_vertex.end());
    //     auto bending_constraint = new BendingConstraint<TinyScalar, TinyConstants>(
    //         vertex_indices, mMesh->mVertices[i], mMesh->getInitialVertices(vertex_indices),  
    //         mParameter.bend);
    //     mConstraints.push_back(bending_constraint);
    //     /* code */
    // }



    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> v(vertices.size()*3);
    for(int i =0;i<vertices.size();i++) {
        v.template block<3,1>(i*3,0) = vertices[i];

    }

    if (this->mIsAttach) {
        for (int i = 0; i < vertices.size(); i++) {
            if (this->mAttachIdx.find(i) != this->mAttachIdx.end()) {
                auto attach_constraint = new AttachmentConstraint<TinyScalar, TinyConstants>(
                    this->mAttachStiffness, i, vertices[i]);
                mConstraints.push_back(attach_constraint);
            }
        }
    }

    for(int i=0;i<this->mMesh->mTriangles.size();i++) {
        this->mContactFace.push_back(this->mMesh->mTriangles[i]);
        for (int j=0;j<3;j++) {
            int idx_f = this->mContactFace[i][j];
            if (std::find(this->mContactIdx.begin(), this->mContactIdx.end(), idx_f) 
                    == this->mContactIdx.end() )
                this->mContactIdx.push_back(idx_f);
        }
    }


    this->m_masses.resize(vertices.size());
    for (std::size_t vId = 0; vId < this->m_masses.size(); ++vId)
    {
        auto triangle_around_vertex_mass_range =
          mCollisionMesh->getTriangleAroundVertexMassRange(vId);
        this->m_masses[vId] = std::accumulate(boost::begin(triangle_around_vertex_mass_range),
                                        boost::end(triangle_around_vertex_mass_range),
                                        TinyConstants::zero()) / 3;

        // std::cout << vId << " " << this->m_masses[vId] << " vMass\n";
    }
    

    world->AddBody(v,mConstraints,this->m_masses,this->mContactIdx,this->mContactFace);
}


template <typename TinyScalar, typename TinyConstants> 
void Dfcloth<TinyScalar, TinyConstants>::
SetMesh(const std::string& path,const Eigen::Transform<TinyScalar, 3, 2>& T)
{
    mScalingMatrix = T;
    this->mMesh = new OBJClothMesh<TinyScalar, TinyConstants>(
        std::string(SOFTCON_DIR)+"/data/"+path,mScalingMatrix
    );


    std::vector<Eigen::Matrix<TinyScalar, 3, 1>> positions;
    std::vector<int> triangles;
    for (int i = 0; i < this->mMesh->mVertices.size(); i++) {
        positions.push_back(this->mMesh->mVertices[i]);
    }
    for (int i = 0; i < this->mMesh->mTriangles.size(); i++) {
        triangles.push_back(this->mMesh->mTriangles[i][0]);
        triangles.push_back(this->mMesh->mTriangles[i][1]);
        triangles.push_back(this->mMesh->mTriangles[i][2]);
        this->mContours.push_back(this->mMesh->mTriangles[i]);
    }
    mCollisionMesh = new CollisionMesh<TinyScalar, TinyConstants>(positions, triangles);

}


template <typename TinyScalar, typename TinyConstants> 
std::vector<std::vector<TinyScalar>> Dfcloth<TinyScalar, TinyConstants>::
getFrictionCoefficientsFromArray(const rapidjson::Value& array)
{
    std::size_t material_number = array.Size();
    std::vector<std::vector<TinyScalar>>  coefficients_table(material_number,
                                                             std::vector<TinyScalar>(material_number));

    for (std::size_t i = 0; i < material_number; ++i)
    {
        if (!array[i].IsArray() || array[i].Size() != material_number)
        {
            throw std::invalid_argument("JSON reader - Error : Friction "
                                        "coefficients must be a square matrix");
        }
        for (std::size_t j = 0; j < material_number; ++j)
        {
            if (!array[i][j].IsDouble())
                throw std::invalid_argument("JSON reader - Error : Friction "
                                            "coefficients must be a square "
                                            "matrix of TinyScalar");
            coefficients_table[i][j] = array[i][j].GetDouble();
        }
    }
    return coefficients_table;
}


// std::vector<std::unique_ptr<Obstacle>>
// getObstaclesFromArray(const rapidjson::Value& array)
// {
//     std::vector<std::unique_ptr<Obstacle>> result;
//     std::transform(array.Begin(), array.End(), std::back_inserter(result), getObstacleFromObject);
//     return result;
// }
#endif