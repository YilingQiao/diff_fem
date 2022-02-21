#ifndef __CONSTRAINT_H__
#define __CONSTRAINT_H__	
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
namespace FEM
{
enum ConstraintType
{
	ATTACHMENT,
	COROTATE,
	LINEAR_MUSCLE,
	BENDING,
	STRAIN,
	TRIANGLE_MUSCLE,
	PNEUMATIC
};
template <typename TinyScalar, typename TinyConstants> 
class Constraint
{
public:
	Constraint(const TinyScalar& stiffness);

	virtual int GetDof() = 0;
	virtual ConstraintType GetType() = 0;

	virtual void	EvaluateJMatrix(int index, std::vector<Eigen::Triplet<TinyScalar>>& J_triplets) =0;
	virtual void	EvaluateLMatrix(std::vector<Eigen::Triplet<TinyScalar>>& L_triplets) =0;
	virtual void 	EvaluateDVector(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& x) = 0;
	virtual void 	GetDVector(int& index,Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1>& d) = 0;
	
	virtual void	fixIndex(int offset) = 0;
// protected:
	TinyScalar mStiffness;

};

template <typename TinyScalar, typename TinyConstants> 
Constraint<TinyScalar, TinyConstants>::
Constraint(const TinyScalar& stiffness)
	:mStiffness(stiffness)
{
	
}
};
#endif