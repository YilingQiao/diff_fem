#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdio.h>
#include <vector>  

#include "../sim/diff/tiny_double_utils.h"
#include "../sim/diff/cppad_utils.h"

#include "../render/SimulationWindow.h"
// #ifdef __APPLE__
// #include <OpenGL/gl.h>
// #include <OpenGL/glu.h>
// #include <GLUT/glut.h>
// #else
// #include <GL/gl.h>
// #include <GL/glu.h>
// #include <GL/glut.h>
// #endif

typedef CppAD::AD<double> Scalar;
typedef CppADUtils<> Utils;

// std::vector<double> mJac;
std::vector<Scalar> mGlobalAx;
bool mDiffState;
bool mDiffAction;
bool mIsBackwardedOnce;
int mTotalStep;
int mDimDiff;
int mDimState;
int mDimAction;
int mCurrentStep;




template <typename TinyScalar, typename TinyConstants>
Environment<TinyScalar, TinyConstants>* init_env(
    int total_step,  
    const std::string& config_path,
    bool do_collision=true,
    bool diff_state=true,
    bool diff_action=false) 
{
    mDiffState = diff_state;
    mDiffAction = diff_action;
    mTotalStep = total_step;
    mCurrentStep = 0;
    mDimDiff = 0; mDimState = 0; mDimAction = 0;
    mIsBackwardedOnce = false;


    Environment<TinyScalar, TinyConstants>* mEnv = 
        new Environment<TinyScalar, TinyConstants>(total_step, config_path, do_collision);
    std::cout << diff_state << diff_action << "!!!!!!!!!!!!!!!\n";
    if (diff_state) {
        mDimState += (mEnv->mSoftWorld->mX.size() + 
            mEnv->mSoftWorld->mV.size()) * mTotalStep;
    }

    mDimDiff = mDimState + mDimAction;

    // for (int i = 0; i < mDimDiff; i++)
    //     mGlobalAx.push_back(TinyConstants::zero());
    // CppAD::Independent(mGlobalAx);
    return mEnv;
}



template <typename TinyScalar, typename TinyConstants>
std::vector<std::vector<double> > forward_act(
    std::vector<double> &q, std::vector<double> &qd, std::vector<double> &jq, 
    std::vector<double> &jqd, std::vector<double> &tau,
    Environment<TinyScalar, TinyConstants>* mEnv,
    bool is_act) 
{
  

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> next_x, next_v, next_jq, next_jqd;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tx(q.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tv(qd.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tjq(jq.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tjqd(jqd.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ttau(tau.size());
    for (int i = 0; i < q.size(); i++) 
        tx[i] = TinyConstants::scalar_from_double(q[i]);
    for (int i = 0; i < qd.size(); i++) 
        tv[i] = TinyConstants::scalar_from_double(qd[i]);
    for (int i = 0; i < jq.size(); i++) 
        tjq[i] = TinyConstants::scalar_from_double(jq[i]);
    for (int i = 0; i < jqd.size(); i++) 
        tjqd[i] = TinyConstants::scalar_from_double(jqd[i]);
    for (int i = 0; i < tau.size(); i++) 
        ttau[i] = TinyConstants::scalar_from_double(tau[i]);

    if (is_act) {
        int cur_phase = mEnv->GetPhase();
        mEnv->SetActions(ttau);
        mEnv->SetPhase(cur_phase+1);
    }

    mEnv->Step(tx, tv, tjq, tjqd, ttau, next_x, next_v, next_jq, next_jqd);


    std::vector<double> ansq, ansqd, ansjq, ansjqd;

    for (int i = 0; i < mEnv->mSoftWorld->mX.size(); i++) {
        ansq.push_back(TinyConstants::getDouble(next_x[i]));
    }
    for (int i = 0; i < mEnv->mSoftWorld->mV.size(); i++) {
        ansqd.push_back(TinyConstants::getDouble(next_v[i]));
    }
    for (int i = 0; i < mEnv->mSoftWorld->mJointQ.size(); i++) {
        ansjq.push_back(TinyConstants::getDouble(next_jq[i]));
    }
    for (int i = 0; i < mEnv->mSoftWorld->mJointQd.size(); i++) {
        ansjqd.push_back(TinyConstants::getDouble(next_jqd[i]));
    }

  return {ansq, ansqd, ansjq, ansjqd};
}     



template <typename TinyScalar, typename TinyConstants>
std::vector<std::vector<double> > backward_act(
    std::vector<double> &q, std::vector<double> &qd, std::vector<double> &jq, 
    std::vector<double> &jqd,
    std::vector<double> &tau, std::vector<double> &dldq, 
    std::vector<double> &dldqd, std::vector<double> &dldjq,
    std::vector<double> &dldjqd,
    Environment<TinyScalar, TinyConstants>* mEnv,
    bool is_act) 
{

    std::vector<TinyScalar> ax(q.size()+qd.size()+jq.size()+jqd.size()+tau.size());
    for (int i = 0; i < ax.size(); i++)
        ax[i] = TinyConstants::zero();

    CppAD::Independent(ax);


    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> next_x, next_v, next_jq, next_jqd;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tx(q.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tv(qd.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tjq(jq.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tjqd(jqd.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ttau(tau.size());
    for (int i = 0; i < q.size(); i++) 
        tx[i] = TinyConstants::scalar_from_double(q[i]) + ax[i];
    for (int i = 0; i < qd.size(); i++) 
        tv[i] = TinyConstants::scalar_from_double(qd[i]) + ax[i+q.size()];
    for (int i = 0; i < jq.size(); i++) 
        tjq[i] = TinyConstants::scalar_from_double(jq[i]) + ax[i+q.size()+qd.size()];
    for (int i = 0; i < jqd.size(); i++) 
        tjqd[i] = TinyConstants::scalar_from_double(jqd[i]) + ax[i+q.size()+qd.size()+jq.size()];
    for (int i = 0; i < tau.size(); i++) 
        ttau[i] = TinyConstants::scalar_from_double(tau[i]) + ax[i+q.size()+qd.size()+jq.size()+jqd.size()];
    

    if (is_act) {
        int cur_phase = mEnv->GetPhase();
        mEnv->SetActions(ttau);
        mEnv->SetPhase(cur_phase+1);
    } else {
        for(int i=0; i<mEnv->mDfobj->mMuscles.size(); i++)
        {
            auto& muscle = mEnv->mDfobj->mMuscles[i];
            for(int j=0;j<muscle->mSegments.size();j++)
                muscle->mSegments[j]->SetActivationLevel(muscle->mActivationLevels[j]);
            }
    }

    mEnv->Step(tx, tv, tjq, tjqd, ttau, next_x, next_v, next_jq, next_jqd);

    
    TinyScalar loss = 0;
    for (int i = 0; i < next_x.size(); i++) 
        loss = loss + TinyConstants::scalar_from_double(dldq[i]) * next_x[i];
    for (int i = 0; i < next_v.size(); i++) 
        loss = loss + TinyConstants::scalar_from_double(dldqd[i]) * next_v[i];
    for (int i = 0; i < next_jq.size(); i++) 
        loss = loss + TinyConstants::scalar_from_double(dldjq[i]) * next_jq[i];
    for (int i = 0; i < next_jqd.size(); i++) 
        loss = loss + TinyConstants::scalar_from_double(dldjqd[i]) * next_jqd[i];

    CppAD::ADFun<double> f(ax, {loss});

    std::vector<double> double_ax(ax.size());
    for (int i = 0; i < ax.size(); i++)
        double_ax[i] = 0.;
    std::vector<double> jac  = f.Jacobian(double_ax); 

    std::vector<double> ansq, ansqd, ansjq, ansjqd, anstau;  
    for (int i = 0; i < next_x.size(); i++) 
        ansq.push_back(jac[i]);
    for (int i = 0; i < next_v.size(); i++) 
        ansqd.push_back(jac[next_x.size() + i]);
    for (int i = 0; i < jq.size(); i++) 
        ansjq.push_back(jac[next_x.size() + next_v.size() + i]);
    for (int i = 0; i < jqd.size(); i++) 
        ansjqd.push_back(jac[next_x.size() + next_v.size() + next_jq.size() + i]);
    for (int i = 0; i < tau.size(); i++) 
        anstau.push_back(jac[next_x.size() + next_v.size() + next_jq.size() + next_jqd.size() + i]);

    return {ansq, ansqd, ansjq, ansjqd, anstau};
}


template <typename TinyScalar, typename TinyConstants>
std::vector<std::vector<double> > forward_step1(
    std::vector<double> &q, std::vector<double> &qd, std::vector<double> &jq, 
    std::vector<double> &jqd, 
    std::vector<double> &tau,
    Environment<TinyScalar, TinyConstants>* mEnv) 
{
  

    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> next_x, next_v, next_jq, next_jqd;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tx(q.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tv(qd.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tjq(jq.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tjqd(jqd.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ttau(tau.size());
    for (int i = 0; i < q.size(); i++) 
        tx[i] = TinyConstants::scalar_from_double(q[i]);
    for (int i = 0; i < qd.size(); i++) 
        tv[i] = TinyConstants::scalar_from_double(qd[i]);
    for (int i = 0; i < jq.size(); i++) 
        tjq[i] = TinyConstants::scalar_from_double(jq[i]);
    for (int i = 0; i < jqd.size(); i++) 
        tjqd[i] = TinyConstants::scalar_from_double(jqd[i]);
    for (int i = 0; i < tau.size(); i++) 
        ttau[i] = TinyConstants::scalar_from_double(tau[i]);
    mEnv->Step(tx, tv, tjq, tjqd, ttau, next_x, next_v, next_jq, next_jqd);

    std::vector<double> ansq, ansqd, ansjq, ansjqd;

    for (int i = 0; i < mEnv->mSoftWorld->mX.size(); i++) {
        ansq.push_back(TinyConstants::getDouble(next_x[i]));
    }
    for (int i = 0; i < mEnv->mSoftWorld->mV.size(); i++) {
        ansqd.push_back(TinyConstants::getDouble(next_v[i]));
    }
    for (int i = 0; i < mEnv->mSoftWorld->mJointQ.size(); i++) {
        ansjq.push_back(TinyConstants::getDouble(next_jq[i]));
    }
    for (int i = 0; i < mEnv->mSoftWorld->mJointQd.size(); i++) {
        ansjqd.push_back(TinyConstants::getDouble(next_jqd[i]));
    }

  return {ansq, ansqd, ansjq, ansjqd};
}     



template <typename TinyScalar, typename TinyConstants>
std::vector<std::vector<double> > backward_step1(
    std::vector<double> &q, std::vector<double> &qd, std::vector<double> &jq, 
    std::vector<double> &jqd, std::vector<double> &tau, std::vector<double> &dldq, 
    std::vector<double> &dldqd, std::vector<double> &dldjq, std::vector<double> &dldjqd,
    Environment<TinyScalar, TinyConstants>* mEnv) 
{

    std::vector<TinyScalar> ax(q.size()+qd.size()+jq.size()+jqd.size()+tau.size());
    for (int i = 0; i < ax.size(); i++)
        ax[i] = TinyConstants::zero();

    CppAD::Independent(ax);


    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> next_x, next_v, next_jq, next_jqd;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tx(q.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tv(qd.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tjq(jq.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tjqd(jqd.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ttau(tau.size());
    for (int i = 0; i < q.size(); i++) 
        tx[i] = TinyConstants::scalar_from_double(q[i]) + ax[i];
    for (int i = 0; i < qd.size(); i++) 
        tv[i] = TinyConstants::scalar_from_double(qd[i]) + ax[i+q.size()];
    for (int i = 0; i < jq.size(); i++) 
        tjq[i] = TinyConstants::scalar_from_double(jq[i]) + ax[i+q.size()+qd.size()];
    for (int i = 0; i < jqd.size(); i++) 
        tjqd[i] = TinyConstants::scalar_from_double(jqd[i]) + ax[i+q.size()+qd.size()+jq.size()];
    for (int i = 0; i < tau.size(); i++) 
        ttau[i] = TinyConstants::scalar_from_double(tau[i]) + ax[i+q.size()+qd.size()+jq.size()+jqd.size()];
    mEnv->Step(tx, tv, tjq, tjqd, ttau, next_x, next_v, next_jq, next_jqd);
    
    TinyScalar loss = 0;
    for (int i = 0; i < next_x.size(); i++) 
        loss = loss + TinyConstants::scalar_from_double(dldq[i]) * next_x[i];
    for (int i = 0; i < next_v.size(); i++) 
        loss = loss + TinyConstants::scalar_from_double(dldqd[i]) * next_v[i];
    for (int i = 0; i < next_jq.size(); i++) 
        loss = loss + TinyConstants::scalar_from_double(dldjq[i]) * next_jq[i];
    for (int i = 0; i < next_jqd.size(); i++) 
        loss = loss + TinyConstants::scalar_from_double(dldjqd[i]) * next_jqd[i];

    CppAD::ADFun<double> f(ax, {loss});

    std::vector<double> double_ax(ax.size());
    for (int i = 0; i < ax.size(); i++)
        double_ax[i] = 0.;
    std::vector<double> jac  = f.Jacobian(double_ax); 

    std::vector<double> ansq, ansqd, ansjq, ansjqd, anstau;  
    for (int i = 0; i < next_x.size(); i++) 
        ansq.push_back(jac[i]);
    for (int i = 0; i < next_v.size(); i++) 
        ansqd.push_back(jac[next_x.size() + i]);
    for (int i = 0; i < jq.size(); i++) 
        ansjq.push_back(jac[next_x.size() + next_v.size() + i]);
    for (int i = 0; i < jqd.size(); i++) 
        ansjqd.push_back(jac[next_x.size() + next_v.size() + next_jq.size() + i]);
    for (int i = 0; i < tau.size(); i++) 
        anstau.push_back(jac[next_x.size() + next_v.size() + next_jq.size() + next_jqd.size() + i]);

    return {ansq, ansqd, ansjq, ansjqd, anstau};
}

template <typename TinyScalar, typename TinyConstants>
std::vector<std::vector<double> > forward_all(
    std::vector<double> &q, std::vector<double> &qd, std::vector<double> &tau,
    Environment<TinyScalar, TinyConstants>* mEnv) 
{

    std::vector<TinyScalar>& ax = mGlobalAx;
    ax.resize((q.size()+qd.size()+tau.size()));
    for (int i = 0; i < ax.size(); i++)
        ax[i] = TinyConstants::zero();

    CppAD::Independent(ax);


    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> next_x, next_v;
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tx(q.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> tv(qd.size());
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ttau(tau.size());
    for (int i = 0; i < q.size(); i++) 
        tx[i] = TinyConstants::scalar_from_double(q[i]) + ax[i];
    for (int i = 0; i < qd.size(); i++) 
        tv[i] = TinyConstants::scalar_from_double(qd[i]) + ax[i+q.size()];
    for (int i = 0; i < tau.size(); i++) 
        ttau[i] = TinyConstants::scalar_from_double(tau[i]) + ax[i+q.size()+qd.size()];

    for (int i = 0; i < mTotalStep; i++) {
        if (i%6==0) {
            int cur_phase = mEnv->GetPhase();
            mEnv->SetActions(ttau);
            mEnv->SetPhase(cur_phase+1);
        }
        mEnv->Step();


        const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x1 
            = mEnv->GetSoftWorld()->GetPositions();
        char filename[50];
        sprintf(filename, "out/%03d.obj",i);
        mEnv->GetDfobj()->OutputSurface(filename, x1);

    }

    next_x = mEnv->mSoftWorld->mX;
    next_v = mEnv->mSoftWorld->mV;
    
    std::vector<double> ansq, ansqd;

    for (int i = 0; i < mEnv->mSoftWorld->mX.size(); i++) {
        ansq.push_back(TinyConstants::getDouble(next_x[i]));
    }
    for (int i = 0; i < mEnv->mSoftWorld->mV.size(); i++) {
        ansqd.push_back(TinyConstants::getDouble(next_v[i]));
    }

    return {ansq, ansqd};
}     



template <typename TinyScalar, typename TinyConstants>
std::vector<std::vector<double> > backward_all(
    const std::vector<double> &q, const std::vector<double> &qd, 
    const std::vector<double> &tau, std::vector<double> &dldq, 
    std::vector<double> &dldqd,
    Environment<TinyScalar, TinyConstants>* mEnv) 
{

    std::vector<TinyScalar>& ax = mGlobalAx;

    TinyScalar loss = 0;
    for (int i = 0; i < q.size(); i++) 
        loss = loss + TinyConstants::scalar_from_double(dldq[i]) * q[i];
    for (int i = 0; i < qd.size(); i++) 
        loss = loss + TinyConstants::scalar_from_double(dldqd[i]) * qd[i];

    CppAD::ADFun<double> f(ax, {loss});

    std::vector<double> double_ax(ax.size());
    for (int i = 0; i < ax.size(); i++)
        double_ax[i] = 0.;
    std::vector<double> jac  = f.Jacobian(double_ax); 

    std::vector<double> ansq, ansqd, anstau;  
    for (int i = 0; i < q.size(); i++) 
        ansq.push_back(jac[i]);
    for (int i = 0; i < qd.size(); i++) 
        ansqd.push_back(jac[q.size() + i]);
    for (int i = 0; i < tau.size(); i++) 
        anstau.push_back(jac[q.size() + qd.size() + i]);
    mEnv->Reset();
    
    return {ansq, ansqd, anstau};
}


template <typename TinyScalar, typename TinyConstants>
std::vector<double> vector2double(
    const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& ev) 
{
    std::vector<double> v;
    for (int i = 0; i < ev.rows(); i++)
        v.push_back(TinyConstants::getDouble(ev[i]));

    return v;
}


template <typename TinyScalar, typename TinyConstants>
Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> double2vector(
    const std::vector<double>& v) 
{
    Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> ev(v.size());
    for (int i = 0; i < v.size(); i++)
        ev[i] = TinyConstants::scalar_from_double(v[i]);

    return ev;
}

template <typename TinyScalar, typename TinyConstants>
void step(Environment<TinyScalar, TinyConstants>* mEnv) {
    static int i = 0;
    mEnv->Step();
    // const Eigen::VectorXd& x1 = mEnvironment->GetSoftWorld()->GetPositions();
    // char filename[50];
    // sprintf(filename, "out/%03d.obj",i);
    // mEnvironment->GetDfobj()->OutputSurface(filename, x1);
    ++i;
}

namespace py = pybind11;

PYBIND11_MODULE(pydifem, m) {
    m.doc() = R"pbdoc(
        differentiable soft physics python plugin
        -----------------------

        .. currentmodule:: pydifem

        .. autosummary::
           :toctree: _generate

    )pbdoc";

    m.def("step", &step<Scalar, Utils>, py::return_value_policy::copy);
        // .def("init_env", ::init_env, py::return_value_policy::copy)
    m.def("init_env", &init_env<Scalar, Utils>, 
        py::arg("total_step"), py::arg("config_path"), py::arg("do_collision")=true, 
        py::arg("diff_state")=true, py::arg("diff_action")=false);
        // .def("run_with_gl", ::run_with_gl, py::return_value_policy::copy)
    ;
        
    // m.def("forward_step", ::forward_step<Scalar, Utils>, py::return_value_policy::copy)
    //     .def("backward_step", ::backward_step<Scalar, Utils>, py::return_value_policy::copy)
    m.def("forward_step1", ::forward_step1<Scalar, Utils>, py::return_value_policy::copy)
        .def("backward_step1", ::backward_step1<Scalar, Utils>, py::return_value_policy::copy)
        .def("forward_all", ::forward_all<Scalar, Utils>, py::return_value_policy::copy)
        .def("backward_all", ::backward_all<Scalar, Utils>, py::return_value_policy::copy)
        .def("forward_act", ::forward_act<Scalar, Utils>, py::return_value_policy::copy)
        .def("backward_act", ::backward_act<Scalar, Utils>, py::return_value_policy::copy);

    m.def("vector2double", ::vector2double<Scalar, Utils>, py::return_value_policy::copy);
    m.def("double2vector", ::double2vector<Scalar, Utils>, py::return_value_policy::copy);

    py::class_<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> >(m, "EigenMatrix")
        .def("rows", &Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::rows)
        .def("cols", &Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::cols);
    py::class_<Eigen::Matrix<Scalar, Eigen::Dynamic, 1> >(m, "EigenVector")
        .def("rows", &Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::rows)
        .def("cols", &Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::cols);

    py::class_<World<Scalar, Utils>,
        std::unique_ptr<World<Scalar, Utils>>>(m, "World")
        .def_readwrite("mX", &World<Scalar, Utils>::mX)
        .def_readwrite("mV", &World<Scalar, Utils>::mV)
        .def_readwrite("mJointQ", &World<Scalar, Utils>::mJointQ)
        .def_readwrite("mJointQd", &World<Scalar, Utils>::mJointQd)
        .def_readwrite("mJointTorque", &World<Scalar, Utils>::mJointTorque)
        .def("GetTimeStep", &World<Scalar, Utils>::GetTimeStep); 

    py::class_<Dfobj<Scalar, Utils>,
        std::unique_ptr<Dfobj<Scalar, Utils>>>(m, "Dfobj")
        .def("SetActivationLevelsAggregate", &Dfobj<Scalar, Utils>::SetActivationLevelsAggregate)
        .def("GetActivationLevelsAggregate", &Dfobj<Scalar, Utils>::GetActivationLevelsAggregate); 

    py::class_<Environment<Scalar, Utils>,
        std::unique_ptr<Environment<Scalar, Utils>>>(m, "Environment")
        .def("SaveObj", &Environment<Scalar, Utils>::SaveObj)
        .def_readwrite("mActions", &Environment<Scalar, Utils>::mActions)
        .def_readwrite("mPhase", &Environment<Scalar, Utils>::mPhase)
        .def_readwrite("mDfobj", &Environment<Scalar, Utils>::mDfobj)
        .def_readwrite("mSoftWorld", &Environment<Scalar, Utils>::mSoftWorld);

    // py::class_<Environment<Scalar,Utils>>(m, "Environment");  
    py::class_<Scalar>(m, "Scalar");
    py::class_<Utils>(m, "Utils")
        .def_static("getDouble", &Utils::getDouble)
        .def_static("scalar_from_double", &Utils::scalar_from_double)
        .def_static("zero", &Utils::zero)
        .def_static("fraction", &Utils::fraction);



  m.attr("__version__") = "dev";
}
