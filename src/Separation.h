#pragma once

#include "EigenTypes.h"

#include <list>

using namespace std;

#ifndef INF
#define INF numeric_limits<double>::infinity()
#endif

class Separation
{
public:
	enum class SeparationEnergy { LOG, QUADRATIC, FLAT_LOG, QUOTIENT, QUOTIENT_NEW };

	Separation();

	void init(int n);
	void value(const MatX2& X, double& f);
	void gradient(const MatX2& X, Vec& g);
	void hessian(const MatX2& X);

	void find_single_hessian(const Vec2& xi, const Vec2& xj, Mat4& h);
	void update_alphas(const Mat& weights, double max_possible);

	SpMat EVvar1, EVvar2, Esep, Esept, V2V, V2Vt;
	SpMat C2C; //Corner to corner
	MatX2 EsepP;

	double Lsep = 1.0, delta = 1.0;
	SeparationEnergy sepEType = SeparationEnergy::QUOTIENT_NEW;

	Vec f_per_pair, f_sep_per_pair;

	// pardiso vectors
	vector<int> II, JJ;
	vector<double> SS;

	// force these uv vertices to be connected more closely, used for gradient
	vector<int> gradient_force_connects;
	// same for function value, to affect the correct index in f_per_row
	// since here its already sorted according to pairs
	vector<int> value_force_connects;

	double force_factor = 10.;

	// weighting indicated by the coloring of the mesh
	// alphas gathered by summing up the factors
	// for each corner force
	Vec connect_alphas;
	// same vars for disconnect
	Vec disconnect_alphas;

	Vec edge_lenghts_per_pair;
	Vec no_seam_constraints_per_pair;
	vector<std::pair<int, int>> pair2ind;
	map<std::pair<int, int>,int> ind2pair;
private:
	Vec EsepP_squared_rowwise_sum;
	Vec EsepP_squared_rowwise_sum_plus_delta;
	

	void flat_log_single_hessian(const Vec2& xi, const Vec2& xj, Mat4& h);
	void make_spd(Mat4& h);
	void add_to_global_hessian(const Mat4& sh, int idx_xi, int idx_xj, int n, list<Tripletd>& htriplets);
	inline int sign(double val);
	inline double dirac(double val);

	void prepare_hessian(int n);
};