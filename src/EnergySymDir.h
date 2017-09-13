#pragma once

#include "EigenTypes.h"
#include "Utils.h"
#include<functional>
#include <Eigen/Core>
#include <Eigen/Sparse>

class DistortionSymDir
{

public:

	/**************************************************************************************************************************/
	//INITIALIZATION 
	DistortionSymDir();

	void init(const MatX3& V, const MatX3i& F, const MatX2& Vs, const MatX3i& Fs);

	void value(const MatX2& X, double& f);

	void gradient(const MatX2& X, Vec& g);

	void hessian(const MatX2& X);

	//loop implementation
	void prepare_hessian(int n);
	/****************************************************************************************************************************/
	double bound=0;
	Eigen::MatrixX3i F;
	Eigen::MatrixX2d V;

	int numV;
	int numE;
	int numS;
	int numF;

	//Jacobian of the parameterization per face
	Eigen::VectorXd a;
	Eigen::VectorXd b;
	Eigen::VectorXd c;
	Eigen::VectorXd d;
	//Eigen::MatrixXd Juv;		//[a,b,c,d]
	//Eigen::MatrixXd invJuv;	//the order and the signs isn't important because the energy is squared anyway thus will be [a,b,c,d]*1/(ad-bc)
	Eigen::VectorXd detJuv;		//(ad-bc)
	Eigen::VectorXd invdetJuv;	//1/(ad-bc)
	Eigen::SparseMatrix<double> DdetJuv_DUV; //jacobian of the function (detJuv) by UV

	//singular values
	Eigen::MatrixX2d s; //Singular values s[0]>s[1]
	Eigen::MatrixX4d v; //Singular vectors 
	Eigen::MatrixX4d u; //Singular vectors 
	Eigen::MatrixXd Dsd[2]; //singular values dense derivatives s[0]>s[1]

	//SVD methods
	bool updateJ(const MatX2& X);
	void UpdateSSVDFunction();
	void ComputeDenseSSVDDerivatives();


	//loop implementation
	inline Eigen::Matrix<double, 6, 6> ComputeFaceConeHessian(const Eigen::Matrix<double,6,1> A1, const Eigen::Matrix<double, 6, 1>& A2, double a1x, double a2x);
	inline Mat6 ComputeConvexConcaveFaceHessian( const Vec6& a1, const Vec6& a2, const Vec6& b1, const Vec6& b2, double aY, double bY, double cY, double dY, const Vec6& dSi, const Vec6& dsi, double gradfS, double gradfs, double HS, double Hs);

	//Energy parts
	//distortion
	Eigen::VectorXd Efi;     //Efi=sum(Ef_dist.^2,2), for data->Efi history

	Eigen::MatrixXi Fuv;                             //F of cut mesh for u and v indices 6XnumF
	Eigen::VectorXd Area;
	Eigen::Matrix3Xd D1d, D2d;						//dense mesh derivative matrices

	Eigen::SparseMatrix<double> a1, a1t, a2, a2t, b1, b1t, b2, b2t;     //constant matrices for cones calcualtion
	Eigen::MatrixXd a1d, a2d, b1d, b2d;					//dense constant matrices for cones calcualtion

//per face hessians vector
   std::vector<Eigen::Matrix<double,6,6>> Hi;
   // pardiso variables
   vector<int> II, JJ;
   vector<double> SS;
};