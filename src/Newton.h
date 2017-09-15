#pragma once

#include "Solver.h"
#include "EigenTypes.h"

#ifdef USE_PARDISO
#include "PardisoSolver.h"
#endif

#include <iostream>
#include <Eigen/SparseCholesky>

using namespace std;

class Newton : public Solver
{
public:
	Newton();

	int step();
	void linesearch();
	bool test_progress();
	void internal_init();
	void internal_update_external_mesh();

private:
	// Wrapper function for flip_avoiding_line_search
	double eval_ls(Mat& x);

	// multiply a std::vector by a constant
	void mult(vector<double>& v, double s);

	// norm of the progress on the mesh
	double diff_norm;

	// Solver that computes Hp = -g
	Eigen::SimplicialLDLT<SpMat> solver; // not used anymore

#ifdef USE_PARDISO
	unique_ptr<PardisoSolver<vector<int>, vector<double>>> pardiso = nullptr;
#else
  	bool needs_init = true;
#endif
	long long int prevTime;
};