#ifndef __WORLD__H__
#define __WORLD__H__
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <vector>
#include <set>
#include <Eigen/Dense>
#include "Constraint/ConstraintHeader.h"
#include "../utils/collisionBrutal.h"
#include "../constants.h"
namespace FEM
{   
// class Constraint;

template <typename TinyScalar, typename TinyConstants> 
class World
{
public:
    World(
        TinyScalar time_step = 1.0/100.0,
        int max_iteration = 30,
        TinyScalar damping_coeff = 0.999,
        bool is_precessing_collision = false
        );
    ~World();
    void                                Initialize();
    void                                Reset();

    void                                AddBody(
                                            const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x0,
                                            const std::vector<Constraint<TinyScalar, TinyConstants>*>& c,
                                            const std::vector<TinyScalar>& masses, 
                                            const std::vector<int>& contactIdx, 
                                            const std::vector<Eigen::Vector3i>& contactFace);

    void                                AddBody(
                                            const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x0,
                                            const std::vector<Constraint<TinyScalar, TinyConstants>*>& c,
                                            const std::vector<TinyScalar>& masses, 
                                            const std::vector<int>& contactIdx, 
                                            const std::vector<Eigen::Vector3i>& contactFace,
                                            const std::vector<std::vector<int>> &rigidBodyIndices, 
                                            const std::vector<int> &nonrigid, 
                                            const std::vector<int> &rigid, 
                                            const std::vector<Link<TinyScalar, TinyConstants>*> &links);

    void                                AddBody(
                                            const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x0,
                                            const std::vector<int>& contactIdx, 
                                            const std::vector<Eigen::Vector3i>& contactFace);

    void                                AddConstraint(Constraint<TinyScalar, TinyConstants>* c);
    void                                RemoveConstraint(Constraint<TinyScalar, TinyConstants>* c);

    void                                TimeStepping(bool isIntegrated = true);
    // void                                TimeSteppingRigid(bool isIntegrated = true);
    void                                UpdatePositionsAndVelocities(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x_n1);

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                        ProjectiveDynamicsMethod();
    void                                PreComputation();

    void                                EvaluateJMatrix(Eigen::SparseMatrix<TinyScalar>& J);
    void                                EvaluateLMatrix(Eigen::SparseMatrix<TinyScalar>& L);
    void                                EvaluateAMatrix();
    void                                EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d);
    void                                FactorizeLDLT(const Eigen::SparseMatrix<TinyScalar>& A,Eigen::SimplicialLDLT<Eigen::SparseMatrix<TinyScalar>>& ldltSolver);

    void                                SetExternalForce(Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> external_force);
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>&             GetExternalForce() {return mExternalForces;};

    const TinyScalar&                       GetTimeStep(){return mTimeStep;};
    const TinyScalar&                       GetTime(){return mTime;};

    const std::vector<Constraint<TinyScalar, TinyConstants>*>&      GetConstraints(){return mConstraints;};
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>&                 GetPositions(){return mX;};
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>&                 GetVelocities(){return mV;};
    const int&                          GetNumVertices(){return mNumVertices;};
    void do_solve_friction(Eigen::Matrix<TinyScalar, 3, 1> &p0, TinyScalar m0, Eigen::Matrix<TinyScalar, 3, 1> &r0, std::vector<Eigen::Matrix<TinyScalar, 3, 1>> &ps, std::vector<TinyScalar> &ms, std::vector<Eigen::Matrix<TinyScalar, 3, 1>> &rs, std::vector<bool> &stick, TinyScalar mu);

    

// private:
    bool                                mIsInitialized;
    bool                                mIsCollision;

    std::vector<Constraint<TinyScalar, TinyConstants>*>         mConstraints;
    int                                 mConstraintDofs;
    int                                 mNumVertices;

    TinyScalar                              mTimeStep;
    TinyScalar                              mTime;
    int                                 mFrame;
    int                                 mMaxIteration, mMaxIteration1;
    TinyScalar                              mDampingCoefficient;

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                        mX,mV, mNonrigidX,mJointQ, mJointQd;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                        mInitX,mInitV;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                        mExternalForces;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                        mQn; // (qt + hvt + h^2M^-1fext)
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                     m_rhs_speed;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                     mJointTorque;

    std::vector<TinyScalar>                 mUnitMass;
    Eigen::SparseMatrix<TinyScalar>         mMassMatrix;
    Eigen::SparseMatrix<TinyScalar>         mInvMassMatrix;
    Eigen::SparseMatrix<TinyScalar>         mIdentityMatrix;

    Eigen::SparseMatrix<TinyScalar>         mJ,mL;
    Eigen::SparseMatrix<TinyScalar>         m_A;
    Eigen::SparseMatrix<TinyScalar>         m_ATA;
    Eigen::SparseMatrix<TinyScalar>         mMh2L;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<TinyScalar>>  
        mDynamicSolver,mQuasiStaticSolver,mDynamicSolverM;

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> mV_without_collision;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> mCollision_V_next;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> mCollision_X_next;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> mPosition_current_next;

    TinyScalar self_collision_tolerance = 0.01; // TODO
    TinyScalar m_self_collision_tol2;

    bool                                m_handle_self_collision = true;
    int m_current_index;
    int m_next_index;

    /// @brief (IO)
    std::vector<TinyScalar> m_rhs_base_times;
    /// @brief (IO)
    std::vector<TinyScalar> m_friction_times;
    /// @brief (IO)
    std::vector<TinyScalar> m_self_friction_times;
    std::vector<TinyScalar> m_self_collision_ordering_times;
    using ForceType = Eigen::Matrix<TinyScalar, 3, -1, Eigen::RowMajor>;
    /// @brief Storage needed to handle self friction
    ForceType m_contact_forces[2];
    /// @brief Storage needed to handle self friction
    ForceType m_self_contact_forces[2];
    /// @brief Storage needed to handle self friction
    ForceType m_self_contact_repercusion_forces[2];
    /// @brief Storage needed to handle self friction
    Eigen::Matrix<TinyScalar, 3, Eigen::Dynamic> m_alpha[2];
    struct SelfCollisionGraphNeighboor
    {
        /**
         * The index of vertex adjacent through this edge.
         */
        std::size_t index;
        /**
         * The index of the self-collision represented by this index.
         */
        std::size_t collision_index;
    };
    /**
     * An adjacency list representing the self-collision graph. The vertices of
     * this graph are the vertices of the mesh. The edges are the self-collision
     * between the vertices.
     * @see SelfCollisionGraphNeighboor
     * @see computeSelfCollisionGraph
     */



    int                                 mNumContactVertices;
    std::map<int, int>                  mMapAll2Contact;
    std::vector<Eigen::Vector3i>        mContactFace, mObsFace;
    std::vector<int>                    mContactIdx;
    CollisionMesh<TinyScalar, TinyConstants>*   mCollisionMesh;
    // Left for children's equations
    // TODO : make it clean
    /// @brief Right hand side of the global step equation
    Eigen::Matrix<TinyScalar, 3, -1, Eigen::RowMajor> m_rhs;
    std::vector<TinyScalar> m_masses;


    CollisionBrutal<TinyScalar, TinyConstants>*  mCollisionDetector;

    void computeCollisionComputationOrder() noexcept;
    std::vector<int> CollisionHandling();
    TinyScalar getVertexMass(size_t vertex_index) const;
    void updateCollisionScene() noexcept;
    void convertCollisionIndexBack();
    void convertCollisionIndex();
    void updateCollisionInIter(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& v, const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x);

    void                                ComputeV0NTr();
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                     ComputeOriX(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &x_n1);
    void                                ComputeSolution(
                                            const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &deltar, 
                                            const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &ori_c, 
                                            Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &hvf, 
                                            Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &s);
    std::vector<std::vector<int>>       mrigidBodies;
    int                                 mDof, mStoreDof;
    TinyScalar                              mAddTorque, mAlpha, mFriction;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic>    QbMmQcLQfm1QcLT, QbMmQcLMfm1QcLT;
    Eigen::SparseMatrix<TinyScalar>         mH2ML, Qf, Qc, Qb, Qc_L, Qb_M, mLf, h2Mf, mC;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                     mcf, mcs, V0;
    std::vector<int>                    mnonrigidIdx, mrigidIdx, mrigidSizes;
    std::vector<Transformation<TinyScalar, TinyConstants> >         T, Tr;
    std::vector<Link<TinyScalar, TinyConstants>*>                  mlinks;


    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>                     v_n1;
    std::vector<Eigen::Matrix<TinyScalar, 3, 1> >        mCollisionMeshPositions, mObsPositions;
    ArcsimMesh                          *mMesh, *mObsMesh;
    Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> mperm;
};




#define EPS 5e-7

template <typename TinyScalar, typename TinyConstants> 
World<TinyScalar, TinyConstants>::
World(
        TinyScalar time_step,
        int max_iteration,
        TinyScalar damping_coeff,
        bool is_precessing_collision)
    :mTimeStep(time_step),
    mFrame(0),
    mMaxIteration(max_iteration),
    mDampingCoefficient(damping_coeff),
    mNumVertices(0),
    mNumContactVertices(0),
    mConstraintDofs(0),
    mIsCollision(is_precessing_collision),
    mIsInitialized(false)
{
    mDof = 0;
    mStoreDof = 0;
}

template <typename TinyScalar, typename TinyConstants> 
World<TinyScalar, TinyConstants>::
~World() {
    delete mCollisionMesh;
    delete mCollisionDetector;
    for(auto c : mConstraints) 
    {
        delete c;
    }
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
Initialize()
{
    // mIsCollision = m_handle_self_collision = false;
    mTime = 0.0;
    mConstraintDofs = 0;

    for(auto c : mConstraints) 
    {
        mConstraintDofs += c->GetDof();
    }

    mV.resize(3*mNumVertices);
    mV.setZero();

    mIdentityMatrix.resize(3*mNumVertices,3*mNumVertices);
    mMassMatrix.resize(3*mNumVertices,3*mNumVertices);
    mInvMassMatrix.resize(3*mNumVertices,3*mNumVertices);
    
    std::vector<Eigen::Triplet<TinyScalar>> i_triplets;
    std::vector<Eigen::Triplet<TinyScalar>> m_triplets;
    std::vector<Eigen::Triplet<TinyScalar>> inv_m_triplets;
    
    i_triplets.reserve(3*mNumVertices);
    m_triplets.reserve(3*mNumVertices);
    inv_m_triplets.reserve(3*mNumVertices);

    for(int i=0;i<mNumVertices;i++)
    {
        m_triplets.push_back(Eigen::Triplet<TinyScalar>(3*i+0,3*i+0,mUnitMass[i]));
        m_triplets.push_back(Eigen::Triplet<TinyScalar>(3*i+1,3*i+1,mUnitMass[i]));
        m_triplets.push_back(Eigen::Triplet<TinyScalar>(3*i+2,3*i+2,mUnitMass[i]));

        inv_m_triplets.push_back(Eigen::Triplet<TinyScalar>(3*i+0,3*i+0,1.0/mUnitMass[i]));
        inv_m_triplets.push_back(Eigen::Triplet<TinyScalar>(3*i+1,3*i+1,1.0/mUnitMass[i]));
        inv_m_triplets.push_back(Eigen::Triplet<TinyScalar>(3*i+2,3*i+2,1.0/mUnitMass[i]));


        i_triplets.push_back(Eigen::Triplet<TinyScalar>(3*i+0,3*i+0,1.0));
        i_triplets.push_back(Eigen::Triplet<TinyScalar>(3*i+1,3*i+1,1.0));
        i_triplets.push_back(Eigen::Triplet<TinyScalar>(3*i+2,3*i+2,1.0));
    }

    mMassMatrix.setFromTriplets(m_triplets.cbegin(), m_triplets.cend());
    mInvMassMatrix.setFromTriplets(inv_m_triplets.cbegin(), inv_m_triplets.cend());
    mIdentityMatrix.setFromTriplets(i_triplets.cbegin(), i_triplets.cend());

    mExternalForces.resize(3*mNumVertices);
    mExternalForces.setZero();

    mQn.resize(3*mNumVertices);
    mQn.setZero();

    mInitX = mX;
    mInitV = mV;

    // collision
    if (m_handle_self_collision)
    {
        m_current_index = 0;
        m_next_index = 1;
        for (int i = 0; i < 2; ++i)
        {
            m_contact_forces[i].setZero(3, mNumContactVertices);
            m_self_contact_forces[i].setZero(3, mNumContactVertices);
            m_self_contact_repercusion_forces[i].setZero(3, mNumContactVertices);
            m_alpha[i].setZero(3, mNumContactVertices);
        } // i
    }     // m_handle_self_collision

    // std::vector<Eigen::Matrix<TinyScalar, 3, 1>> positions;
    std::vector<int> triangles;

    // for (int i = 0; i < mContactIdx.size(); i++) {
    //  int id = mContactIdx[i];
    //  Eigen::Matrix<TinyScalar, 3, 1> v(mX[3*id], mX[3*id+1], mX[3*id+2]);
    //  positions.push_back(v);
    // }
    for (int i = 0; i < mContactFace.size(); i++) {
        triangles.push_back(mContactFace[i][0]);
        triangles.push_back(mContactFace[i][1]);
        triangles.push_back(mContactFace[i][2]);
    }

    // mCollisionMesh = new CollisionMesh<TinyScalar, TinyConstants>(positions, triangles);
    // mCollisionDetector = new CollisionBrutal<TinyScalar, TinyConstants>(mCollisionMesh, true);   
    mCollisionMesh = new CollisionMesh<TinyScalar, TinyConstants>(mCollisionMeshPositions, triangles);
    mObsMesh = new ArcsimMesh();
    mObsMesh->Initialize(
        helper::to_eigen_double<TinyScalar, TinyConstants>(mObsPositions), 
        helper::to_eigen_double<TinyScalar, TinyConstants>(mObsPositions), 
        mObsFace);
    mMesh = new ArcsimMesh();
    mMesh->Initialize(
        helper::to_eigen_double<TinyScalar, TinyConstants>(mCollisionMeshPositions), 
        helper::to_eigen_double<TinyScalar, TinyConstants>(mCollisionMeshPositions), 
        mContactFace);
    mCollisionDetector = new CollisionBrutal<TinyScalar, TinyConstants>(mCollisionMesh, true); 
    mCollisionDetector->mObsMesh = mObsMesh;  
    mCollisionDetector->mMesh = mMesh; 

    mCollisionDetector->m_generalized_positions.resize(mNumContactVertices*3);
    mCollisionDetector->m_generalized_speeds.resize(mNumContactVertices*3);
    mCollisionDetector->mCollisionMesh->m_positions.resize(mNumContactVertices*3);
    mCollisionDetector->mCollisionMesh->m_speed.resize(mNumContactVertices*3);
    m_rhs = Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(
        3, mContactIdx.size());
    
    mCollisionDetector->m_handle_self_collision = m_handle_self_collision;

    mCollision_V_next.resize(mNumContactVertices*3);
    mCollision_X_next.resize(mNumContactVertices*3);
    std::cout<<"m_rhs rows: "<<m_rhs.rows()<<std::endl;
    std::cout<<"m_rhs cols: "<<m_rhs.cols()<<std::endl;
    // collision
    // rigid bodies
    T.resize(mrigidBodies.size());
    Tr.resize(mrigidBodies.size());
    mJointQ.resize(mStoreDof);
    mJointQd.resize(mStoreDof);
    for (int i = 0; i < mrigidBodies.size(); ++i)
        mlinks[i]->WriteT(mJointQ, mJointQd);

    mJointTorque.resize(std::max(0, int(mrigidBodies.size()-1) ));
    mJointTorque.setZero();

    // rigid bodies

    // PreComputation();
    PreComputation();

    mIsInitialized = true;

    std::cout<<"Total degree of freedom : "<<mX.rows()<<std::endl;
    std::cout<<"Total constraints : "<<mConstraints.size()<<std::endl;
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
Reset()
{
    mX = mInitX;
    mV = mInitV;
    mTime = 0.0;
}

template <typename TinyScalar, typename TinyConstants> 
TinyScalar World<TinyScalar, TinyConstants>::
getVertexMass(size_t vertex_index) const
{
    return mH2ML.coeff(mContactIdx[vertex_index]*3,mContactIdx[vertex_index]*3);
    // return mUnitMass[mContactIdx[vertex_index]];
  // return 1.0;
  // return mCollisionDetector->mCollisionMesh->m_masses[vertex_index];
}


// add cloth body
template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
AddBody(
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x0,
    const std::vector<Constraint<TinyScalar, TinyConstants>*>& c,
    const std::vector<TinyScalar>& masses, 
    const std::vector<int>& contactIdx, 
    const std::vector<Eigen::Vector3i>& contactFace)
{
    std::cout << contactFace.size() << " addbody contactFace\n";
    std::cout << contactIdx.size() << " addbody contactIdx\n";
    
    for (const auto& cidx : contactIdx) {
        int id = cidx + mNumVertices;
        if (std::find(mContactIdx.begin(), mContactIdx.end(), id) 
                == mContactIdx.end() ) {
            mContactIdx.push_back(id);
            mMapAll2Contact.insert(std::make_pair(id, mNumContactVertices++));
            mCollisionMeshPositions.push_back(x0.template segment<3>(cidx*3));
        }
    }

    for (const auto& f : contactFace) {
        int f0, f1, f2;
        f0 = mMapAll2Contact[f[0]+mNumVertices];
        f1 = mMapAll2Contact[f[1]+mNumVertices];
        f2 = mMapAll2Contact[f[2]+mNumVertices];
        mContactFace.push_back(Eigen::Vector3i(f0, f1, f2));
    }
    
    int nv=(x0.rows()/3);
    for (int i = 0; i < nv*3; ++i)
        mnonrigidIdx.push_back(mNumVertices*3 + i);
    mNumVertices+=nv;

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tmpX = mX;
    // tmpX.resize(mNumVertices*3);
    mX.resize(mNumVertices*3);

    mX.head(tmpX.rows()) = tmpX;
    mX.tail(x0.rows()) = x0;

    for (auto con : c) {
        con->fixIndex(mNumVertices - nv);
    }
    mConstraints.insert(mConstraints.end(),c.begin(),c.end());

    for(int i=0;i<nv;i++){
        mUnitMass.push_back(masses[i]);
    }

    if(mIsInitialized)
        Initialize();
}

//add obstacle body
template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
AddBody(
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x0,
    const std::vector<int>& contactIdx, 
    const std::vector<Eigen::Vector3i>& contactFace)
{
    int offset = mObsPositions.size();
    for (const auto& cidx : contactIdx) {
        // mContactIdx.push_back(-1);
        // mCollisionMeshPositions.push_back(x0.segment<3>(cidx*3));
        mObsPositions.push_back(x0.template segment<3>(cidx*3));
    }
    // mNumContactVertices += contactIdx.size();

    for (const auto& f : contactFace) {
        int f0, f1, f2;
        f0 = f[0]+offset;
        f1 = f[1]+offset;
        f2 = f[2]+offset;
        mObsFace.push_back(Eigen::Vector3i(f0, f1, f2));
    }

    if(mIsInitialized)
        Initialize();
}


//add rtq8 soft body
template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
AddBody(
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x0,
    const std::vector<Constraint<TinyScalar, TinyConstants>*>& c,
    const std::vector<TinyScalar>& masses, 
    const std::vector<int>& contactIdx, 
    const std::vector<Eigen::Vector3i>& contactFace,
    const std::vector<std::vector<int>> &rigidBodyIndices, 
    const std::vector<int> &nonrigid, 
    const std::vector<int> &rigid, 
    const std::vector<Link<TinyScalar, TinyConstants>*> &links)
{
    printf("AddBody rigid bodies\n");
    int offset = mNumVertices, rigidoffset = mrigidIdx.size()/3;
    for (int i = 0; i < rigidBodyIndices.size(); ++i) {
        auto tmp = rigidBodyIndices[i];
        for (int i = 0; i < tmp.size(); ++i)
            tmp[i] += offset;
        Link<TinyScalar, TinyConstants> *link = links[i];
        link->bodyIdx = mrigidBodies.size();
        mlinks.push_back(link);

        link->rigidoffset = rigidoffset;
        rigidoffset += tmp.size();

        link->dofIdx = mDof;
        link->storedofIdx = mStoreDof;
        if (link->parent == NULL) {
            mDof += 6;
            mStoreDof += 12;
        }
        else {
            mDof += 1;
            mStoreDof += 1;
        }

        mrigidBodies.push_back(tmp);
        mrigidSizes.push_back(tmp.size());
        link->indices = tmp;
    }
    for (int idx : nonrigid)
        mnonrigidIdx.push_back(offset * 3 + idx);
    for (int idx : rigid)
        mrigidIdx.push_back(offset * 3 + idx);
    int k = 0;
    for (auto &idxs : mrigidBodies) {
        for (int i = 0; i < idxs.size(); ++i) {
            if (mrigidIdx[k] != idxs[i]*3 || mrigidIdx[k+1] != idxs[i]*3+1 || mrigidIdx[k+2] != idxs[i]*3+2)
                std::cout << "???" << mrigidIdx[k] << " " << idxs[i]<< std::endl;
            k += 3;
        }
    }

    // collision
    std::cout << contactFace.size() << " addbody contactFace\n";
    std::cout << contactIdx.size() << " addbody contactIdx\n";
    
    for (const auto& cidx : contactIdx) {
        int id = cidx + mNumVertices;
        if (std::find(mContactIdx.begin(), mContactIdx.end(), id) 
                == mContactIdx.end() ) {
            mContactIdx.push_back(id);
            mMapAll2Contact.insert(std::make_pair(id, mNumContactVertices++));
            mCollisionMeshPositions.push_back(x0.template segment<3>(cidx*3));
        }
    }

    for (const auto& f : contactFace) {
        int f0, f1, f2;
        f0 = mMapAll2Contact[f[0]+mNumVertices];
        f1 = mMapAll2Contact[f[1]+mNumVertices];
        f2 = mMapAll2Contact[f[2]+mNumVertices];
        mContactFace.push_back(Eigen::Vector3i(f0, f1, f2));
    }
    // collision end
    int nv=(x0.rows()/3);
    mNumVertices+=nv;

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tmpX = mX;
    // tmpX.resize(mNumVertices*3);
    mX.resize(mNumVertices*3);

    mX.head(tmpX.rows()) = tmpX;
    mX.tail(x0.rows()) = x0;

    for (auto con : c) {
        con->fixIndex(mNumVertices - nv);
    }
    mConstraints.insert(mConstraints.end(),c.begin(),c.end());

    for(int i=0;i<nv;i++){
        mUnitMass.push_back(masses[i]);
    }

    if(mIsInitialized)
        Initialize();


}





template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
AddConstraint(Constraint<TinyScalar, TinyConstants>* c)
{
    mConstraints.push_back(c);
    if(mIsInitialized){
        mConstraintDofs = 0;
        for(auto c : mConstraints){
            mConstraintDofs += c->GetDof();
        }
        PreComputation();
        // PreComputation();
    }
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
RemoveConstraint(Constraint<TinyScalar, TinyConstants>* c)
{
    bool isRemoved = false;
    for(int i=0;i<mConstraints.size();i++)
    {
        if(mConstraints[i]==c) {
            mConstraints.erase(mConstraints.begin()+i);
            isRemoved = true;
            break;
        }
    }

    if(isRemoved) {
        if(mIsInitialized) {
            mConstraintDofs = 0;
            for(auto c : mConstraints){
                mConstraintDofs += c->GetDof();
            }
            PreComputation();
        }
    }
}

// template <typename TinyScalar, typename TinyConstants> 
// void World<TinyScalar, TinyConstants>::
// TimeStepping(bool isIntegrated)
// {
//  if(!mIsInitialized) {
//      std::cout<<"Engine not initialized."<<std::endl;
//      return;
//  }

//     Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> x_n1(mNumVertices*3);
//     Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> v_n1(mNumVertices*3);

//     mV_without_collision = mV + mTimeStep*(mInvMassMatrix*mExternalForces);

//  // (qt + hvt + h^2M^-1fext)
//  mQn = mX + mTimeStep*mV_without_collision;
//  mPosition_current_next = mQn;

//     auto tmp = mTimeStep*(mInvMassMatrix*mExternalForces);


//     x_n1=ProjectiveDynamicsMethod();
//     // v_n1=ProjectiveDynamicsMethod();

//     // UpdatePositionsAndVelocities(v_n1);
//     UpdatePositionsAndVelocities(x_n1);
//  // mV *= mDampingCoefficient;
    
//     updateCollisionScene();
//  if(isIntegrated)
//  {   
//      mTime += mTimeStep;
//      mFrame++;
//  }   
// }

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
TimeStepping(bool isIntegrated)
{
    if(!mIsInitialized) {
        std::cout<<"Engine not initialized."<<std::endl;
        return;
    }
    for (int i = 0; i < mrigidBodies.size(); ++i)
        mlinks[i]->ReadT(mJointQ, mJointQd);

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> x_n1(mNumVertices*3);

    // (qt + hvt + h^2M^-1fext)
    // mExternalForces[0] = 1;
    // std::cout << mInvMassMatrix << " mInvMassMatrix\n";
    // std::cout << mV[1] << " mV\n";
    // std::cout << mExternalForces[1] << " mExternalForces\n";

    mV_without_collision = mV + mTimeStep*(mInvMassMatrix*mExternalForces);

    // std::cout << mV[1] << " mV\n";
    // std::cout << mExternalForces[1] << " mExternalForces\n";
    // std::cout << mV_without_collision[1] << " mV_without_collision\n";

    // (qt + hvt + h^2M^-1fext)
    mQn = mX + mTimeStep*mV_without_collision;
    updateCollisionScene();
    // std::cout << mQn[1] << " mQn\n";
    x_n1 = ProjectiveDynamicsMethod();
    // std::cout << x_n1[1] << " x_n1\n";
    // std::cout << mX[1] << " mX\n";

    UpdatePositionsAndVelocities(x_n1);

    mV *= mDampingCoefficient;
    // std::cout << mV[1] << " mV 2\n";
    updateCollisionScene();

    for (int i = 0; i < mrigidBodies.size(); ++i)
        mlinks[i]->WriteT(mJointQ, mJointQd);

    if(isIntegrated)
    {   
        mTime += mTimeStep;
        mFrame++;
    }   
}



template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
updateCollisionScene() noexcept
{
    for (int i = 0; i < mContactIdx.size(); i++) {
        int id = mContactIdx[i];
        for (int j = 0; j < 3; j++) {
            mCollisionDetector->m_generalized_positions[3*i+j] = mX[3*id+j];
            mCollisionDetector->m_generalized_speeds[3*i+j] = mV[3*id+j];
            mCollisionDetector->mCollisionMesh->m_positions[3*i+j] = mX[3*id+j];
            mCollisionDetector->mCollisionMesh->m_speed[3*i+j] = mV[3*id+j];
        }
        mCollisionDetector->mMesh->nodes[i]->x0 = Vec3(
            TinyConstants::getDouble(mX[3*id+0]) ,
            TinyConstants::getDouble(mX[3*id+1]) ,
            TinyConstants::getDouble(mX[3*id+2]));
    }
   
    // Cleaning the forces for the next step
    if (m_handle_self_collision)
    {
        m_current_index = 0u;
        m_next_index = 1u;
        for (unsigned int i = 0u; i < 2u; ++i)
        {
            m_contact_forces[i].setZero();
            m_self_contact_forces[i].setZero();
            m_self_contact_repercusion_forces[i].setZero();
            m_alpha[i].setZero();
        }
    }
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
UpdatePositionsAndVelocities(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x_n1)
{
    // mV = x_n1;
    // mX = mX + mV * mTimeStep;

    mV = (x_n1 - mX) * (1.0 / mTimeStep);
    mX = x_n1;
}

// template <typename TinyScalar, typename TinyConstants> 
// Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> World<TinyScalar, TinyConstants>::
// ProjectiveDynamicsMethod()
// {
//  Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> x_n1(3*mNumVertices);
//     Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> x_n1_new(3*mNumVertices);
//     Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> x_n1_new_v(3*mNumVertices);
//  Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> b(3*mNumVertices);
//  Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> d(3*mConstraintDofs);
//     Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> rhs_speed(3*mNumVertices);
//     Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> v_n1(3*mNumVertices);
//     Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> v_n1_new(3*mNumVertices);
//  d.setZero();
//  b= (1.0/(mTimeStep*mTimeStep))*mMassMatrix*mQn;

//     TinyScalar h2 = mTimeStep*mTimeStep;

//  x_n1 = mQn;
//     v_n1 = mV_without_collision;

//  if(mIsCollision) {
//      mCollisionDetector->BaseUpdateCollisions(mCollision_X_next);
//         updateCollisionInIter(v_n1, x_n1);
//  }

//  int i;
//     TinyScalar err = 100.;
//  for(i=0; i<mMaxIteration; i++) {


//      EvaluateDVector(x_n1,d);

//         m_rhs_speed = mMassMatrix * mV_without_collision +
//             mTimeStep * (mJ * d - mL * mX);

//         if(mIsCollision) {
//             convertCollisionIndex();
//             CollisionHandling();
//             convertCollisionIndexBack();
//         }

//         v_n1_new = mDynamicSolver.solve(m_rhs_speed/h2);
//         x_n1_new = mX + v_n1_new*mTimeStep;

//         if(mIsCollision) {
//             updateCollisionInIter(v_n1_new, x_n1_new);
//         }
  
//         err = (x_n1_new - x_n1).norm()/x_n1.size();
//      if(err < EPS) {
//          break;
//      }   
//      x_n1 = x_n1_new;
//  }
//     // if(mIsCollision) {
//     //     updateCollisionInIter(v_n1_new, x_n1_new);
//     // }

//     std::cout << "error " << err << " ";
//     if(err < EPS)
//         std::cout << "good!\n";
//     else
//         std::cout << "not converge!\n";

//  return x_n1;
// }

template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> World<TinyScalar, TinyConstants>::
ProjectiveDynamicsMethod()
{
    static int num_max = 0, num_iter = 0;
    static double total_time = 0;

    int numUnknown = mnonrigidIdx.size() + mDof;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> x_n1(numUnknown);
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ori_x(3*mNumVertices); 
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ori_x_old(3*mNumVertices);
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> x_n1_new(numUnknown);
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> b(3*mNumVertices);
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> d(3*mConstraintDofs);
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ori_c, s, xf, xf_old(mnonrigidIdx.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> hvf(mnonrigidIdx.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> deltar(mnonrigidIdx.size());
    deltar.setZero();

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> rhs_speed(3*mNumVertices);
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> v_n1(3*mNumVertices);
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> v_n1_new(3*mNumVertices);
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> prev_mrhsspeed(3*mNumVertices);

    mNonrigidX.resize(mnonrigidIdx.size());

    d.setZero();
    prev_mrhsspeed.setZero();
    b = (1.0 / (mTimeStep * mTimeStep)) * mMassMatrix * mQn;
    // std::cout << mTimeStep << " mTimeStep\n";
    // std::cout << mQn[1] << " mQn\n";
    // std::cout << b[1] << " b\n";
    TinyScalar h2 = mTimeStep*mTimeStep;
    v_n1 = mV_without_collision;


    x_n1 = (mperm * mQn).template segment(0, mnonrigidIdx.size());
    xf_old = x_n1;
    mNonrigidX = (mperm * mX).template segment(0,mnonrigidIdx.size());

    if(mIsCollision) {
        // mCollisionDetector->BaseUpdateCollisions(mCollision_X_next);
        updateCollisionInIter(mV_without_collision, mQn);
    }
            // std::cout << "\tmv0 "<<((mTimeStep*mV_without_collision).template segment<3>(353*3)).transpose() << std::endl;

    ori_x_old.fill(0.);
    for (int i = 0; i < mrigidBodies.size(); ++i) 
        mlinks[i]->guessNext();
    TinyScalar err;
    for(int i = 0; i < mMaxIteration; i++) {
        ComputeV0NTr();
        ori_x = ComputeOriX(x_n1);

        // std::cout << "ori_x "<< i << " " << ori_x[1] << std::endl;

        TinyScalar err = (ori_x - ori_x_old).norm()/ori_x.size();
        if (err < EPS) {
            // std::cout << "iter i0=" << i << std::endl;
            break;
        }

        ori_x_old = ori_x;
        EvaluateDVector(ori_x,d);


        ori_c = b + mJ * d; 
        // std::cout << "b "<< i << " " << b[1] << std::endl;
        // std::cout << "d "<< i << " " << d[1] << std::endl;

        // std::cout << "hvf0 "<< i << " " << hvf[1] << std::endl;
        // std::cout << "deltar "<< i << " " << deltar[1] << std::endl;
        // std::cout << "ori_c "<< i << " " << ori_c[1] << std::endl;
        ComputeSolution(deltar, ori_c, hvf, s);
            // std::cout << "\thvf0 "<<(hvf.template segment<3>(353*3)).transpose() << std::endl;

        xf = mNonrigidX + hvf;
        // std::cout << "mNonrigidX "<< i << " " << mNonrigidX[1] << std::endl;
        // std::cout << "hvf "<< i << " " << hvf[1] << std::endl;
        x_n1_new << xf, s;
        x_n1 = x_n1_new;
    }
    ComputeV0NTr();
    // std::cout << x_n1_new[1] << " x_n1_new\n";
    ori_x = ComputeOriX(x_n1_new);

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> old_hvf = hvf;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> old_oric = ori_c;

    // std::cout << mIsCollision << " mIsCollision\n";
    if (mIsCollision) {
        mCollisionDetector->ClearCollisions();
        ori_x_old.setZero();
        while (true) {
            // printf("mIsCollision\n");
            auto tmpv = (ori_x - mX) / mTimeStep;
            updateCollisionInIter(tmpv, ori_x);
            for (int i = 0; i < mContactIdx.size(); i++) {
                mCollisionDetector->mMesh->nodes[i]->x = Vec3(
                    TinyConstants::getDouble( mCollision_X_next[3*i+0]),
                    TinyConstants::getDouble( mCollision_X_next[3*i+1]),
                    TinyConstants::getDouble( mCollision_X_next[3*i+2]));
            }
            if (mCollisionDetector->BaseUpdateCollisions(mCollision_X_next))
                break;
            // std::cout << "collision start!" << std::endl;
            ++num_max;
            num_iter += 100;
            bool flag = true;
            for (int i = 0; i < mMaxIteration1; ++i) {
                ComputeV0NTr();
                ori_x = ComputeOriX(x_n1);
                err = (ori_x - ori_x_old).cwiseAbs().maxCoeff();//.norm()/ori_x.size();//
                // std::cout << "err1 " << err << std::endl;
                // std::cout << "\tori_x "<<(ori_x.template segment<3>(mContactIdx[86]*3)).transpose() << std::endl;
                if (err < EPS && flag) {
                    --num_max;
                    num_iter += i - 100;
                    // std::cout << "iter i1=" << i << std::endl;
                    break;
                }
                EvaluateDVector(ori_x,d);
                auto tmpv = (ori_x - mX) / mTimeStep;
                updateCollisionInIter(tmpv, ori_x);
                // std::cout << (tmpv.template segment<3>(mContactIdx[38]*3)).transpose() << std::endl;
                // std::cout << ori_x.segment<3>(mContactIdx[689]*3).transpose() << std::endl;

                // m_rhs_speed = mMassMatrix * mV_without_collision;
                m_rhs_speed = mMassMatrix * mV_without_collision +
                    mTimeStep * (mJ * d - mL * mX);
                auto tmp = mTimeStep * mC * tmpv * mTimeStep;
                // std::cout <<"m_rhs_speed bef "<< (m_rhs_speed.template segment<3>(mContactIdx[656]*3)).transpose() << std::endl;
                m_rhs_speed -= tmp;
                // std::cout <<"m_rhs_speed aft "<< (m_rhs_speed.template segment<3>(mContactIdx[656]*3)).transpose() << std::endl;
                // std::cout <<"ori_x "<< (ori_x.template segment<3>(mContactIdx[656]*3)).transpose() << std::endl;
                if (i == 0)
                    prev_mrhsspeed = m_rhs_speed;
                // m_rhs_speed = mperm * m_rhs_speed;
                // m_rhs_speed.segment(0, mnonrigidIdx.size()) -= tmp;
                // m_rhs_speed = mperm.transpose() * m_rhs_speed;
                // for (int i = 0; i < mnonrigidIdx.size(); i++) {
                //     m_rhs_speed[mnonrigidIdx[i]] -= tmp[i];
                // }
                convertCollisionIndex();
                std::vector<int> cindices = CollisionHandling();
                flag = true;
                // for (int i : cindices) {
                //     auto tmp = (ori_x - ori_x_old).template segment<3>(mContactIdx[i]*3);
                //     // std::cout << i << " " << mContactIdx[i]*3 << " "<< tmp.transpose() << std::endl;
                //     if (tmp.norm() > 1e-5)
                //         flag = false;
                // }
                // std::cout << flag << std::endl;
                // m_rhs_speed.setZero();
                convertCollisionIndexBack();
                //mrhsspeed is delta
                // auto diff = m_rhs_speed - prev_mrhsspeed;
                // std::cout << diff.norm() << std::endl;
                // std::cout << diff.segment<3>(mContactIdx[406]*3).transpose() << std::endl;
                // int idx = 0;
                // for (int i = 0; i < mNumVertices*3; ++i)
                //     if (std::abs(diff[i])>std::abs(diff[idx]))
                //         idx = i;
                // std::cout << "argmax= "<<idx<<" : "<<diff[idx]<<" "<< m_rhs_speed[idx]<<" "<<prev_mrhsspeed[idx]<< std::endl;
                TinyScalar alpha = mAlpha;
                prev_mrhsspeed = (m_rhs_speed * alpha + prev_mrhsspeed * (1 - alpha));
                ori_c = b + mJ * d + prev_mrhsspeed / mTimeStep;
                // std::cout << prev_mrhsspeed.norm() << std::endl;
                // std::cout << (old_oric - ori_c).norm() << std::endl;
                // std::cout <<"m_rhs_speed "<< (m_rhs_speed.template segment<3>(mContactIdx[656]*3)).transpose() << std::endl;
                // ComputeSolution_collision(ori_c, hvf, hvf, s);
                // deltar = mTimeStep * mLf * (hvf - old_hvf);
                // old_hvf = hvf;

            
                // ComputeSolution(deltar, ori_c, hvf, s);
                hvf = mDynamicSolverM.solve((mperm * prev_mrhsspeed / mTimeStep).segment(0, mnonrigidIdx.size()));
                // std::cout <<"hvf "<< (hvf.template segment<3>(mContactIdx[449]*3)).transpose() << std::endl;
                // std::cout <<"hvf "<< (hvf.template segment<3>(mContactIdx[515]*3)).transpose() << std::endl;
                // std::cout <<"hvf "<< (hvf.template segment<3>(mContactIdx[517]*3)).transpose() << std::endl;
                // std::cout <<"hvf "<< (hvf.template segment<3>(mContactIdx[656]*3)).transpose() << std::endl;
                // std::cout <<"mass "<< mH2ML.coeff(mContactIdx[656]*3,mContactIdx[656]*3) << std::endl;
                xf = mNonrigidX + hvf;
                x_n1_new << xf, s;
                x_n1 = x_n1_new;
                ori_x_old = ori_x;
            }
        }     
    }


    // ComputeV0NTr();
    // ori_x = ComputeOriX(x_n1_new);
    return ori_x;
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
convertCollisionIndex() {
    for (int i = 0; i < mContactIdx.size(); i++) {
        int id = mContactIdx[i];
        m_rhs.template block<3,1>(0,i) = m_rhs_speed.template block<3,1>(id*3,0);
    }
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
updateCollisionInIter(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& v, const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x) {
    for (int i = 0; i < mContactIdx.size(); i++) {
        int id = mContactIdx[i];
        mCollision_V_next.template block<3,1>(i*3,0) = v.template block<3,1>(id*3,0);
        mCollision_X_next.template block<3,1>(i*3,0) = x.template block<3,1>(id*3,0);
    }
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
convertCollisionIndexBack() {
    for (int i = 0; i < mContactIdx.size(); i++) {
        int id = mContactIdx[i];
        m_rhs_speed.template block<3,1>(id*3,0) = m_rhs.template block<3,1>(0,i);
    }
}

// template <typename TinyScalar, typename TinyConstants> 
// void World<TinyScalar, TinyConstants>::
// PreComputation()
// {
//  EvaluateJMatrix(mJ);
//  EvaluateLMatrix(mL);
//     EvaluateAMatrix();


//  Eigen::SparseMatrix<TinyScalar> H2ML = (1.0/(mTimeStep*mTimeStep))*mMassMatrix+mL;
//     printf("\n############### mL\n");
//     for (int i = 0; i < 10; i++) { // for debug testing
//         printf("(%f) ", mL.coeff(3*i, 3*i));
//     }
//     printf("\n############### mMassMatrix\n");
//     for (int i = 0; i < 10; i++) { // for debug testing
//         printf("(%f) ", mMassMatrix.coeff(3*i, 3*i));
//     }
//     printf("\n############### H2ML\n");
//     Eigen::SparseMatrix<TinyScalar> h2 = H2ML*(mTimeStep*mTimeStep);
//     for (int i = 0; i < 10; i++) { // for debug testing
//         printf("(%f) ", h2.coeff(3*i, 3*i));
//     }

//     mMh2L = h2;
//  FactorizeLDLT(H2ML,mDynamicSolver);
//  FactorizeLDLT(mL,mQuasiStaticSolver);
// }


template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
PreComputation()
{
    EvaluateJMatrix(mJ);
    EvaluateLMatrix(mL);
    EvaluateAMatrix();
    // mH2ML = (1.0/(mTimeStep*mTimeStep))*mMassMatrix+mL;
    h2Mf = (1.0/(mTimeStep*mTimeStep))*mMassMatrix;
    mH2ML = h2Mf + mL;

    //split
    Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(3*mNumVertices);
    Eigen::SparseMatrix<TinyScalar> tmp = mH2ML;
    auto lmd0 = [](const Eigen::Index& row, const Eigen::Index& col, const auto&){return row==col;};
    auto lmd1 = [](const Eigen::Index& row, const Eigen::Index& col, const auto&){return row!=col;};
    
    const int nnon = mnonrigidIdx.size(), nrig = mrigidIdx.size();
    Eigen::VectorXi tmpp(3*mNumVertices);
    for (int i = 0; i < nnon; ++i)
        tmpp[i] = mnonrigidIdx[i];
    for (int i = 0; i < nrig; ++i)
        tmpp[i + nnon] = mrigidIdx[i];
    //column-major p
    for (int i = 0; i < 3*mNumVertices; ++i)
        perm.indices()[tmpp[i]] = i;
    tmp = mH2ML.twistedBy(perm);
    Qf = tmp.block(0,0,nnon, nnon);
    Qc_L = tmp.block(nnon, 0,nrig, nnon);
    Qb_M = tmp.block(nnon, nnon,nrig, nrig);

    tmp = mH2ML;
    tmp.prune(lmd1);
    mC = tmp*1;

    mMh2L = mH2ML*(mTimeStep*mTimeStep);
    FactorizeLDLT(Qf,mDynamicSolver);

    tmp = mH2ML;
    tmp.prune(lmd0);
    FactorizeLDLT(tmp,mDynamicSolverM);


        

    QbMmQcLQfm1QcLT = Qb_M - Qc_L * mDynamicSolver.solve(
        Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic>(Qc_L.transpose()));
    // FactorizeLDLT(H2ML,mDynamicSolver);
    // FactorizeLDLT(mL,mQuasiStaticSolver);
    mperm = perm;
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
EvaluateJMatrix(Eigen::SparseMatrix<TinyScalar>& J) 
{
    J.resize(3*mNumVertices,3*mConstraintDofs);
    std::vector<Eigen::Triplet<TinyScalar>> J_triplets;
    int index = 0;
    for(int i =0;i<mConstraints.size();i++)
    {
        mConstraints[i]->EvaluateJMatrix(index,J_triplets);
        index+=mConstraints[i]->GetDof();
    }
    J.setFromTriplets(J_triplets.cbegin(), J_triplets.cend());
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
EvaluateAMatrix() 
{
    Eigen::SparseMatrix<TinyScalar> tmpB = mL;
    m_ATA.resize(mNumContactVertices, mNumContactVertices);
    // m_A.resize(mConstraintDofs, mNumContactVertices);


    for (int i = 0; i < mNumContactVertices; i++) {
        for (int j = 0; j < mNumContactVertices; j++) {
            int id = mContactIdx[i];
            int jd = mContactIdx[j];
            if (id == jd) continue;
            // if (id == -1 || jd == -1) continue;
            if (tmpB.coeff(3*id, 3*jd) == 0)
                continue;
            m_ATA.coeffRef(i,j) = tmpB.coeff(3*id, 3*jd);
        }
    }


}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
EvaluateLMatrix(Eigen::SparseMatrix<TinyScalar>& L) 
{
    L.resize(3*mNumVertices,3*mNumVertices);

    std::vector<Eigen::Triplet<TinyScalar>> L_triplets;

    for(auto c : mConstraints)
        c->EvaluateLMatrix(L_triplets);

    L.setFromTriplets(L_triplets.cbegin(), L_triplets.cend());
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d) 
{
    d.resize(mConstraintDofs*3);
    int n = mConstraints.size();
// #pragma omp parallel for
    for(int i=0;i<n;i++)
    {
        mConstraints[i]->EvaluateDVector(x);
    }
    int index = 0;
    for(auto& c : mConstraints)
    {
        c->GetDVector(index,d);
    }
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
FactorizeLDLT(const Eigen::SparseMatrix<TinyScalar>& A,Eigen::SimplicialLDLT<Eigen::SparseMatrix<TinyScalar>>& ldltSolver)
{
    Eigen::SparseMatrix<TinyScalar> A_prime = A;
    ldltSolver.analyzePattern(A_prime);
    ldltSolver.factorize(A_prime);
    TinyScalar reg = 1E-6;
    bool success = true;
    while (ldltSolver.info() != Eigen::Success)
    {
        reg *= 10;
        A_prime = A + reg*mIdentityMatrix;
        ldltSolver.factorize(A_prime);
        success = false;
    }
    if (!success)
        std::cout << "factorize failure (damping : " << reg<<" )"<<std::endl;
    // else
    //     std::cout << "factorize success\n";
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
SetExternalForce(Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> external_force)
{
    mExternalForces = external_force;
}

//  collision ----------------------------------------------------------------------------------------------------------

#define BLOCK(m, id) ((m).col((id)))
#define LVAL(m, id, cmp) ((m)((cmp), (id)))

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
do_solve_friction(
    Eigen::Matrix<TinyScalar, 3, 1> &p0, 
    TinyScalar m0, 
    Eigen::Matrix<TinyScalar, 3, 1> &r0, 
    std::vector<Eigen::Matrix<TinyScalar, 3, 1>> &ps, 
    std::vector<TinyScalar> &ms, 
    std::vector<Eigen::Matrix<TinyScalar, 3, 1>> &rs, 
    std::vector<bool> &stick, 
    TinyScalar mu) {
    //compute stick p and m
    Eigen::Matrix<TinyScalar, 3, 1> pstick = p0;
    TinyScalar mstick = m0;
    int n = stick.size();
    for (int k = 0; k < n; ++k)
        if (stick[k]) {
            pstick += ps[k];
            mstick += ms[k];
        }
    //compute slide p and m and fn
    Eigen::Matrix<TinyScalar, 3, 1> pslide(0.,0.,0.);
    TinyScalar mslide = 0, fn = 0;
    for (int k = 0; k < n; ++k)
        if (!stick[k]) {
            pslide += ps[k];
            mslide += ms[k];
            fn += rs[k][0];
        }
    //compute stick force:
    Eigen::Matrix<TinyScalar, 3, 1> v = pstick / mstick;
    r0[1] = v[1] * m0 - p0[1];
    r0[2] = v[2] * m0 - p0[2];
    for (int k = 0; k < n; ++k)
        if (stick[k]) {
            rs[k][1] = v[1] * ms[k] - ps[k][1];
            rs[k][2] = v[2] * ms[k] - ps[k][2];
        }
    if (mslide > 0) {
        //compute relative v direction, and normalize
        Eigen::Matrix<TinyScalar, 3, 1> rel_v = pslide / mslide - pstick / mstick;
        TinyScalar t_len = sqrt(rel_v[1]*rel_v[1] + rel_v[2]*rel_v[2]);
        rel_v /= t_len;
        //apply sliding frictional force to stuck verts
        TinyScalar mult = mu*fn/(mslide+mstick);
        r0[1] += rel_v[1] * mult * m0;
        r0[2] += rel_v[2] * mult * m0;
        for (int k = 0; k < n; ++k)
            if (stick[k]) {
                rs[k][1] += rel_v[1] * mult * ms[k];
                rs[k][2] += rel_v[2] * mult * ms[k];
            }
        //apply sliding frictional force to slide verts
        Eigen::Matrix<TinyScalar, 3, 1> abs_v0 = (p0 + r0) / m0;
        for (int k = 0; k < n; ++k)
            if (!stick[k]) {
                Eigen::Matrix<TinyScalar, 3, 1> rel_v = ps[k] / ms[k] - abs_v0;
                TinyScalar t_len = sqrt(rel_v[1]*rel_v[1] + rel_v[2]*rel_v[2]);
                rel_v /= t_len;
                TinyScalar mult = -mu*rs[k][0] * ms[k] / (ms[k] + mstick);
                rs[k][1] = mult*rel_v[1];
                rs[k][2] = mult*rel_v[2];
            }
    }
    //sanity check: sum of frictional forces should be one
    Eigen::Matrix<TinyScalar, 3, 1> rsum = r0;
    for (int k = 0; k < n; ++k)
        rsum += rs[k];
    assert(rsum.norm() < 1e-8);
}

template <typename TinyScalar, typename TinyConstants> 
std::vector<int> World<TinyScalar, TinyConstants>::
CollisionHandling() {
    const unsigned int nbCollisions = mCollisionDetector->m_collisions_infos.size();
    const unsigned int nbSelfCollisions = mCollisionDetector->m_self_collisions_infos.size();
    // Abort if there is no collisions
    
    // std::cout << nbSelfCollisions << " nbSelfCollisions \n";
    // std::cout << nbCollisions << " nbCollisions \n";

    if ((nbCollisions == 0) && (!m_handle_self_collision || (nbSelfCollisions == 0)))
    {
        // m_rhs.setZero();
        return {};
    }
    std::vector<int> ans;
    TIMER_START(friction);


    // Reconstruct the forces without any friction force

    Eigen::Matrix<TinyScalar, 3, Eigen::Dynamic> forces = m_rhs;
    // m_rhs.setZero();
    // const auto v = Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic>::Map(mCollision_V_next.data(), 3, mNumContactVertices);

        // std::cout << mNumContactVertices << " mNumContactVertices\n";
        // std::cout << mCollision_V_next.rows() << "!" << mCollision_V_next.cols() << "\n";
        // std::cout << v.rows() << " 1 " << v.cols() << "\n";
        // std::cout << m_ATA.rows() << " 2 " << m_ATA.cols() << "\n";
        // std::cout << forces.rows() << " 3 " << forces.cols() << "\n";

        // std::cout << v.block<3,1>(0,0) << " v\n";
        // std::cout << mCollision_V_next.block<3,1>(0,0) << " v\n";
        // exit(0);
// #pragma omp parallel for
    // for (size_t cmp = 0u; cmp < 3u; ++cmp)
    // {
    //     forces.row(cmp) -=
    //       // std::pow(m_time_step, 2) * (getATA() * v.row(cmp).transpose()).transpose();
    //       //// ! better perf !
    //       std::pow(mTimeStep, 2) * v.row(cmp) * m_ATA;
    // }


    // Estimated friction force
// #pragma omp parallel for
    for (unsigned int cId = 0u; cId < nbCollisions; ++cId)
    {
        // Data
        const CollisionInfo<TinyScalar, TinyConstants> & collision_info = 
            mCollisionDetector->m_collisions_infos[cId];
        const size_t vId = collision_info.vertex_index;
        const Eigen::Matrix<TinyScalar, 3, 3>& local_basis = mCollisionDetector->m_local_contact_basis[cId];
        const TinyScalar mu = collision_info.friction_coefficient;


        // Converting the rhs in the local basis
        // Contact force also sees the previous self collision forces
        Eigen::Matrix<TinyScalar, 3, 1> forceLoc;
        if (!m_handle_self_collision)
        {
            forceLoc = local_basis.transpose() *
                       (BLOCK(forces, vId));// - getVertexMass(vId) * collision_info.speed);
        }
        else
        {
            forceLoc = local_basis.transpose() *
                       (BLOCK(forces, vId) +
                        BLOCK(m_self_contact_forces[m_current_index], vId)
                        );
        }

        // Estimating the reaction
        ans.push_back(vId);
        if (forceLoc[0] > 0.)
        {
            // take-off
            continue;
        }
        Eigen::Matrix<TinyScalar, 3, 1> r(0., 0., 0.);
        // rN must prevent the penetration
        r[0] = -forceLoc[0];

        // rT try to prevent the sliding
        // Sticking
        r[1] = -forceLoc[1];
        r[2] = -forceLoc[2];

        const TinyScalar rT_norm = sqrt(r[1] * r[1] + r[2] * r[2]);
        if (rT_norm > mu * r[0])
        {
            // but gets stuck at the border of the cone
            // Sliding
            r[1] *= mu * r[0] / rT_norm;
            r[2] *= mu * r[0] / rT_norm;
        }
        
        // Converting r in the global frame
        r = local_basis * r;

        if ((!m_handle_self_collision) || (nbSelfCollisions == 0))
        {
            // Adding it directly to the rhs
            for (unsigned int cmp = 0u; cmp < 3u; ++cmp)
            {
                //#pragma omp atomic update Not needed : max 1 contact per vertex
                LVAL(m_rhs, vId, cmp) = r[cmp];
            } // cmp
        }

        if (m_handle_self_collision)
        {
            // Storing it
            for (unsigned int cmp = 0u; cmp < 3u; ++cmp)
            {
                LVAL(m_contact_forces[m_next_index], vId, cmp) = r[cmp];
            } // cmp
        }     // self

    } // cId


    const TinyScalar duration_friction = TIMER_DURATION(friction, microseconds);
    m_friction_times.push_back(duration_friction);
#ifdef TIMER_PRINT
    std::cout << "# Rhs friction : " << duration_friction << " µs" << std::endl;
#endif // TIMER_PRINT


    if ((m_handle_self_collision) && (nbSelfCollisions > 0))
    {
        TIMER_START(self_friction);
        // Estimated self friction force
        // for (size_t level = 0u; level < mCollisionDetector->m_collision_computation_order.size(); ++level)
        {
            // const std::vector<size_t>& collisions_level_ids = mCollisionDetector->m_collision_computation_order[level];

            std::vector<SelfForceToAdd<TinyScalar, TinyConstants> > forces_to_add;
            // forces_to_add.resize(collisions_level_ids.size());
            forces_to_add.resize(nbSelfCollisions);

// #pragma omp parallel for
            for (unsigned int i = 0u; i < forces_to_add.size(); ++i)
            // for (unsigned int i = 0u; i < collisions_level_ids.size(); ++i)
            {
                // const size_t scId = collisions_level_ids[i];
                const size_t scId = i;
                // Data
                const SelfCollisionInfo<TinyScalar, TinyConstants>& self_collision_info = 
                    mCollisionDetector->m_self_collisions_infos[scId];
                const size_t vId = self_collision_info.vertex_index;
                std::vector<size_t> vId_list_old = self_collision_info.vertex_index_list;
                const std::array<int, 3> fId = self_collision_info.face_indices;
                const Eigen::Matrix<TinyScalar, 3, 3>& local_basis = mCollisionDetector->m_local_self_contact_basis[scId];
                const Eigen::Matrix<TinyScalar, 3, 1>& alpha = self_collision_info.barycentric_coordinates;
                TinyScalar mu = mFriction; //self_collision_info.friction_coefficient;

                //compute p and m and transfer to local basis
                TinyScalar m0 = getVertexMass(fId[0]) + getVertexMass(fId[1]) + getVertexMass(fId[2]);
                Eigen::Matrix<TinyScalar, 3, 1> f0(0.,0.,0.), r0(0.,0.,0.);
                for (int i = 0; i < 3; ++i) {
                    f0 += (BLOCK(forces, fId[i]) + BLOCK(m_contact_forces[m_next_index], fId[i]));
                }
                f0 = local_basis.transpose() * f0;
                std::vector<size_t> vId_list;
                for (int k : vId_list_old) {
                    Eigen::Matrix<TinyScalar, 3, 1> f = local_basis.transpose() * (BLOCK(forces, k) + BLOCK(m_contact_forces[m_next_index], k));
                    Eigen::Matrix<TinyScalar, 3, 1> rel_v =  f / getVertexMass(k) - f0 / m0;
                    if (rel_v[0] <= 0)
                        vId_list.push_back(k);
                }
                std::vector<TinyScalar> ms;
                std::vector<Eigen::Matrix<TinyScalar, 3, 1>> fs, rs(vId_list.size());
                TinyScalar smv = 0., sm = 0.;
                for (int k : vId_list) {
                    ms.push_back(getVertexMass(k));
                    sm += ms.back();
                    fs.push_back(local_basis.transpose() * (BLOCK(forces, k) + BLOCK(m_contact_forces[m_next_index], k)));
                    smv += fs.back()[0];
                }
                //compute fn
                for (int k = 0; k < vId_list.size(); ++k) {
                    TinyScalar fn = (f0[0] + smv) * ms[k] / (m0 + sm) - fs[k][0];
                    // Eigen::Matrix<TinyScalar, 3, 1>
                    rs[k].fill(0);
                    rs[k][0] += fn;
                    r0[0] -= fn;
                }
                //assume all stick from the beginning and remove invalid ones
                std::vector<bool> stick;
                for (int k = 0; k < vId_list.size(); ++k)
                    stick.push_back(true);
                while (1) {
                    //output r0, rs
                    do_solve_friction(f0, m0, r0, fs, ms, rs, stick, mu);
                    bool changed = false;
                    for (int k = 0; k < vId_list.size(); ++k)
                        if (stick[k] && rs[k][1]*rs[k][1]+rs[k][2]*rs[k][2] > mu*rs[k][0]*rs[k][0]) {
                            changed = true;
                            stick[k] = false;
                        }
                    if (!changed)
                        break;
                    // std::cout << "changed!" << std::endl;
                }
                //transfer back to world frame
                r0 = local_basis * r0;
                Eigen::Matrix<TinyScalar, 3, 1> v0 = (local_basis* f0 + r0) / m0;
                // std::cout << "expected " << v0.transpose() << std::endl;
                for (int i = 0; i < 3; ++i) {
                    m_self_contact_forces[m_next_index].col(fId[i]) += v0 * getVertexMass(fId[i]) - (BLOCK(forces, fId[i]) + BLOCK(m_contact_forces[m_next_index], fId[i]));
                }
                for (int k = 0; k < vId_list.size(); ++k) {
                    rs[k] = local_basis * rs[k];
                    m_self_contact_forces[m_next_index].col(vId_list[k]) += rs[k];
                }
                
            
            } // scId


// #pragma omp parallel for
//             for (size_t i = 0u; i < collisions_level_ids.size(); ++i)
//             {
//                 const size_t scId = collisions_level_ids[i];
//                 const Eigen::Matrix<TinyScalar, 3, 1> prev_force = mCollisionDetector->m_remember_self_contact_forces.col(scId);
//                 const SelfForceToAdd<TinyScalar, TinyConstants>& new_force = forces_to_add[i];

//                 for (size_t cmp = 0u; cmp < 3u; ++cmp)
//                 {
//                     // Remove from current
// // #pragma omp atomic update
//                     LVAL(m_self_contact_forces[m_current_index], new_force.id_plus, cmp) -=
//                       prev_force[cmp];
// // #pragma omp atomic update
//                       for (int k = 0; k < 3; ++k)
//                     LVAL(m_self_contact_forces[m_current_index], new_force.id_minus[k], cmp) +=
//                       prev_force[cmp] * new_force.alpha[k];
//                     // LVAL(m_self_contact_forces[m_current_index], new_force.id_minus, cmp) +=
//                       // prev_force[cmp];

// // #pragma omp atomic update
//                     LVAL(m_self_contact_forces[m_next_index], new_force.id_plus, cmp) +=
//                       new_force.force[cmp];
// // #pragma omp atomic update
//                       for (int k = 0; k < 3; ++k)
//                     LVAL(m_self_contact_forces[m_next_index], new_force.id_minus[k], cmp) -=
//                       new_force.force[cmp] * new_force.alpha[k];
//                     // LVAL(m_self_contact_forces[m_next_index], new_force.id_minus, cmp) -=
//                     //   new_force.force[cmp];

// // #pragma omp atomic write
//                     mCollisionDetector->m_remember_self_contact_forces(cmp, scId) = new_force.force[cmp];
//                 } // cmp
//             }     // scId

        } // level

        // Add all the forces to the RHS
        // std::cout << "???" << m_self_contact_forces[m_next_index].col(86).transpose()<<std::endl;
        m_rhs += m_contact_forces[m_next_index] + m_self_contact_forces[m_next_index]
          //+ m_self_contact_repercusion_forces[m_next_index]
          ;

        const TinyScalar duration_self_friction = TIMER_DURATION(self_friction, microseconds);
        m_self_friction_times.push_back(duration_self_friction);
#ifdef TIMER_PRINT
        std::cout << "# Rhs self friction : " << duration_self_friction << " µs" << std::endl;
#endif // TIMER_PRINT
    }

    // Move on next step
    if (m_handle_self_collision)
    {
        m_current_index = m_next_index;
        m_next_index = (m_next_index) ? 0u : 1u;
        m_contact_forces[m_next_index].setZero();
        m_self_contact_forces[m_next_index].setZero();
        // m_self_contact_repercusion_forces[m_next_index].setZero();
        m_alpha[m_next_index].setZero();
    }
    return ans;
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
ComputeSolution(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &deltar, 
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &ori_c, 
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &hvf, 
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &s)
{
    // TIMER_START(total)
    //tildaV
    int curr = 0, curc = 0;
    Eigen::SparseMatrix<TinyScalar> tildaV(mrigidIdx.size(), mDof);
    std::vector<Eigen::Triplet<TinyScalar>> trips;
    for (int i = 0; i < mrigidBodies.size(); ++i)
        if (mlinks[i]->parent == NULL)
            mlinks[i]->ComputeTildaV(trips, V0, mInitX);
    tildaV.setFromTriplets(trips.begin(), trips.end());


    Qc = tildaV.transpose() * Qc_L;

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tmp1(mnonrigidIdx.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tmp2(mrigidIdx.size());
    auto tmp = mperm * ori_c;
    tmp1 = tmp.template segment(0, mnonrigidIdx.size());
    tmp2 = tmp.template segment(mnonrigidIdx.size(), mrigidIdx.size());
    // std::cout << "tmp2 size " << tmp2.size() << "\n";


    mcf = tmp1 - Qc_L.transpose() * V0 - Qf*mNonrigidX;
    mcs = tildaV.transpose() * (tmp2 - Qb_M * V0 - Qc_L * mNonrigidX);
    // std::cout << mJointTorque[0] << " " 
    //     << mJointTorque[1] << " "
    //     << mJointTorque[2] << " mJointTorque1\n";
    // std::cout << "mJointTorque size " << mJointTorque.size() << "\n";
    // std::cout << mcs[6] << " " 
    //     << mcs[7] << " "
    //     << mcs[8] << " mcs\n";
    for (int i = 0; i < mJointTorque.size(); ++i) {
        mcs[6+i] += mJointTorque[i];
    }
    // std::cout << mcs[6] << " " 
    //     << mcs[7] << " "
    //     << mcs[8] << " mcs2\n";
    
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> S = tildaV.transpose() * QbMmQcLQfm1QcLT * tildaV;
    // Eigen::Matrix<TinyScalar, Eigen::Dynamic, Eigen::Dynamic> Sinv = S.inverse();
    // TIMER_START(hvf)


        
    if (mrigidIdx.size() > 0) {
        Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> solve_result = mDynamicSolver.solve(-mcf);
    
        Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tmp = mcs + Qc * solve_result;
  
        s = S.fullPivLu().solve(tmp);    
      
        hvf = mDynamicSolver.solve(mcf - Qc.transpose() * s);
    } else {
        hvf = mDynamicSolver.solve(mcf);
        Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ttt = mcf;
        ttt.setZero();
        ttt = mDynamicSolver.solve(ttt);


    }
 
    for (int i = 0; i < mrigidBodies.size(); ++i)
        mlinks[i]->UpdateT(s);
}

template <typename TinyScalar, typename TinyConstants> 
Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> World<TinyScalar, TinyConstants>::
ComputeOriX(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &x_n1)
{
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ans;
    ans.resize(3*mNumVertices);
    for (int i = 0; i < mnonrigidIdx.size(); ++i)
        ans[mnonrigidIdx[i]] = x_n1[i];
    for (int i = 0; i < mrigidIdx.size(); ++i)
        ans[mrigidIdx[i]] = V0[i];
    return ans;
}

template <typename TinyScalar, typename TinyConstants> 
void World<TinyScalar, TinyConstants>::
ComputeV0NTr()
{
    V0.resize(mrigidIdx.size());
    // // std::cout << T.size() << std::endl;
    // for (int i = 0; i < T.size(); ++i) {
    //  Eigen::JacobiSVD<Eigen::Matrix<TinyScalar, 3, 3>> svd(T[i].template block<3,3>(0,0), Eigen::ComputeFullU | Eigen::ComputeFullV);
    //  Tr[i].template block<3,3>(0,0) = svd.matrixU()*svd.matrixV().transpose();
    //  Tr[i].template block<3,1>(0,3) = T[i].template block<3,1>(0,3);
    // // std::cout <<  Tr[i] << std::endl;
    // }
    // for (int i = 0, j = 0; i < mrigidBodies.size(); ++i) {
    //  auto body = mrigidBodies[i];
    //  Eigen::Matrix34d &cur_Tr = Tr[i];
    //  for (int k = 0; k < body.size(); ++k) {
    //      Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> V = mInitX.segment<3>(body[k]*3), vh(4);
    //      vh<<V,1;
    //      Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> v0 = cur_Tr * vh;
    //      V0.segment<3>(j) = v0;
    //      j += 3;
    //  }
    // }
    for (int i = 0; i < mrigidBodies.size(); ++i) 
        if (mlinks[i]->parent == NULL)
            mlinks[i]->ComputeV0NTr(V0, mInitX);
}


#undef EPS

};
#endif