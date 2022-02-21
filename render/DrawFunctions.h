#ifndef __DRAW_FUNCTIONS_H__
#define __DRAW_FUNCTIONS_H__
// #include "../sim/Dfskeleton.h"
#include "../sim/Dfobj.h"
#include "../sim/Dfobstacle.h"
#include "../sim/fem/World.h"
#include "../sim/Muscle.h"
#include "../sim/diff/tiny_double_utils.h"
#include "DrawPrimitives.h"
namespace GUI
{
	void DrawWorld();
    void DrawObs(Dfobstacle<double, DoubleUtils>* Dfobs,const Eigen::Vector3d& eye);
    void DrawCharacter(Dfobj<double, DoubleUtils>* Dfobj,const Eigen::VectorXd& x,const Eigen::Vector3d& eye);
    void DrawCollisionMesh(World<double, DoubleUtils>* world,const Eigen::Vector3d& eye);
	void DrawMuscles(const std::vector<Muscle<double, DoubleUtils>*>& muscles,const Eigen::VectorXd& x);
	void DrawActivations(const double& x, const double& y,const double& length,const double& height,
        Muscle<double, DoubleUtils>* muscle1,Muscle<double, DoubleUtils>* muscle2);
	void DrawOBJ();
};
#endif