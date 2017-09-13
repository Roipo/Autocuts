#include "Separation.h"
#include "autodiff.h"

#include <chrono>
#include <thread>
#include <iostream>
#include <igl/cat.h>

using namespace std;

Separation::Separation()
{
}

void Separation::init(int n)
{
	// wedges energy
	igl::cat(1, EVvar1, EVvar2, Esep);
	Esept = Esep.transpose();
// 	Triplets C2CTriplets;
	for (int i = 0; i < Esept.outerSize(); ++i)
	{
		// no inner loop because there are only 2 nnz values per col
		SpMat::InnerIterator it(Esept, i);
		int idx_xi = it.row();
		int idx_xj = (++it).row();
		pair<int, int> p(idx_xi, idx_xj);
		pair2ind.push_back(p);
		ind2pair.emplace(p, i);
		ind2pair.emplace(std::make_pair(p.second,p.first) , i);
	}

	V2Vt = V2V.transpose();
	C2C = V2Vt*V2V;
 	SpMat I(V2V.cols(), V2V.cols());
 	I.setIdentity();
	C2C -= I;
	C2C.prune(0,0);
	prepare_hessian(n);

	connect_alphas = Vec::Zero(Esep.rows());
	disconnect_alphas = Vec::Ones(Esep.rows());

	edge_lenghts_per_pair = Vec::Ones(Esept.cols());
	no_seam_constraints_per_pair = Vec::Zero(Esept.cols());
}

void Separation::value(const MatX2& X, double& f)
{
	EsepP = Esep * X;
	EsepP_squared_rowwise_sum = EsepP.array().pow(2.0).rowwise().sum();
	EsepP_squared_rowwise_sum_plus_delta = EsepP_squared_rowwise_sum.array() + delta;
	switch (sepEType)
	{
		case SeparationEnergy::LOG:			
			//f = (EsepP_squared_rowwise_sum_plus_delta.array().log() * Lsep).sum();
			f_per_pair = (EsepP_squared_rowwise_sum_plus_delta.array().log() * Lsep);
			break;
		case SeparationEnergy::QUADRATIC:
			//f = EsepP_squared_rowwise_sum.sum();
			f_per_pair = EsepP_squared_rowwise_sum;
			break;
		case SeparationEnergy::FLAT_LOG:
		{ // Verified: f = sum over Lsep * log((||xi-xj||^2 / (||xi-xj||^2 + delta)) + 1) (1 xi-xj pair per col (1 & -1) in Esept
			// Store per-pair value for finding the maximal value in current setup
			f_per_pair = Lsep * (EsepP_squared_rowwise_sum.cwiseQuotient(EsepP_squared_rowwise_sum_plus_delta).array() + 1.0).log();
			break;
		}
		case SeparationEnergy::QUOTIENT:
		case SeparationEnergy::QUOTIENT_NEW:
			f_per_pair = EsepP_squared_rowwise_sum.cwiseQuotient(EsepP_squared_rowwise_sum_plus_delta);
			break;
		default:
			assert(false && "Unimplemented separation energy");
	}
	// store values before taking painting into account
	f_sep_per_pair = f_per_pair;

	// add attraction force from painting
	// alpha * ||xi - xj||^2
	f_per_pair += (connect_alphas + no_seam_constraints_per_pair).cwiseProduct(EsepP_squared_rowwise_sum);

	// apply distraction force from painting
	// f -> alpha * f
	f_per_pair = f_per_pair.cwiseProduct(disconnect_alphas);

	// add edge length factor
	f_per_pair = f_per_pair.cwiseProduct(edge_lenghts_per_pair);

	// if a pair shall not be a seam, it should have a high value
	//f_per_pair = f_per_pair.cwiseProduct(no_seam_constraints_per_pair);

	// sum everything up
	f = f_per_pair.sum();
}

void Separation::gradient(const MatX2& X, Vec& g)
{
	MatX2 ge;
	switch (sepEType)
	{
		case SeparationEnergy::LOG:
			ge = 2.0 * Esep.transpose() * Lsep * EsepP_squared_rowwise_sum_plus_delta.cwiseInverse().asDiagonal() * EsepP;
			break;
		case SeparationEnergy::QUADRATIC:
			ge = 2.0 * Esep.transpose() * EsepP;
			break;
		case SeparationEnergy::FLAT_LOG:
		{	// Verified: dxi = (2 * delta * (xi - xj)) / ((2*||xi-xj||^2 + delta) * (||xi-xj||^2 + delta)) for xi and dxj = -dxi
			Vec d_vec = Vec::Constant(EsepP_squared_rowwise_sum.rows(), delta);
			Vec two_x_plus_a = (2.0 * EsepP_squared_rowwise_sum) + d_vec;
			Vec d = d_vec.cwiseQuotient(EsepP_squared_rowwise_sum_plus_delta.cwiseProduct(two_x_plus_a));
			ge = 2.0 * Esep.transpose() * d.asDiagonal() * EsepP;
			break;
		}
		case SeparationEnergy::QUOTIENT:
		case SeparationEnergy::QUOTIENT_NEW:
		{
			Vec d_vec = Vec::Constant(EsepP_squared_rowwise_sum.rows(), delta);
			Vec x_plus_d = EsepP_squared_rowwise_sum + d_vec;
			Vec d = d_vec.cwiseQuotient(x_plus_d.cwiseAbs2());
			Vec dconn_e_disconn = (d + connect_alphas + no_seam_constraints_per_pair).cwiseProduct(edge_lenghts_per_pair).cwiseProduct(disconnect_alphas);
			ge = 2.0 * Esept * dconn_e_disconn.asDiagonal() * EsepP;
			break;
		}
		default:
			assert(false && "Unimplemented separation energy");
	}
	g = Eigen::Map<Vec>(ge.data(), 2.0 * ge.rows(), 1);
}

void Separation::hessian(const MatX2& X)
{
	int n = X.rows();
	int threads = omp_get_max_threads();
#pragma omp parallel for num_threads(threads)
	for (int i = 0; i < Esept.outerSize(); ++i)
	{ // no inner loop because there are only 2 nnz values per col
		int tid = omp_get_thread_num();
		Vec2 xi, xj;
		Mat4 sh;
		int idx_xi, idx_xj, factor;
		SpMat::InnerIterator it(Esept, i);
		idx_xi = it.row();
		factor = it.value();
		idx_xj = (++it).row();
		xi = X.row(idx_xi);
		xj = X.row(idx_xj);
		find_single_hessian(xi, xj, sh);
		sh *= factor;
		// add the additional factors like coloring and edge splitting/merging
		Mat4 Esep4;
		Esep4 << 1, 0, -1, 0,
			0, 1, 0, -1,
			-1, 0, 1, 0,
			0, -1, 0, 1;
		sh += Esep4 * (connect_alphas(i) + no_seam_constraints_per_pair(i));
		sh *= edge_lenghts_per_pair(i);
		sh *= disconnect_alphas(i);
		//sh *= no_seam_constraints_per_pair(i);

		int ind = 10 * i;
		for (int a = 0; a < 4; ++a)
		{
			for (int b = 0; b <= a; ++b)
			{
				SS[ind++] = sh(b, a);
			}
		}
	}
}

void Separation::prepare_hessian(int n)
{
	II.clear();
	JJ.clear();
	auto PushPair = [&](int i, int j) { II.push_back(i); JJ.push_back(j); };
	for (int i = 0; i < Esept.outerSize(); ++i)
	{
		SpMat::InnerIterator it(Esept, i);
		int idx_xi = it.row();
		int idx_xj = (++it).row();
		// The indices in the small hessians are setup like this:
		// xi, xi+n, xj, xj+n from top to bottom and left to right
		// we traverse only the upper diagonal of each 4x4 hessian
		// and thus store 10 values, gathered in column order.
		// First column
		PushPair(idx_xi,			idx_xi);
		// Second column
		PushPair(idx_xi,			idx_xi + n);
		PushPair(idx_xi + n,	idx_xi + n);
		// Third column
		PushPair(idx_xi,			idx_xj);
		//PushPair(idx_xi + n,	idx_xj);
		PushPair(idx_xj, idx_xi + n);
		PushPair(idx_xj,			idx_xj);
		// Fourth column
		PushPair(idx_xi,			idx_xj + n);
		PushPair(idx_xi + n,	idx_xj + n);
		PushPair(idx_xj,			idx_xj + n);
		PushPair(idx_xj + n,	idx_xj + n);
	}
	SS = vector<double>(II.size(), 0.);
}

void Separation::find_single_hessian(const Vec2& xi, const Vec2& xj, Mat4& h)
{
	bool speedup = true;
	Vec2 dx = xi - xj;
	Vec4 dxx;
	dxx << dx, -dx;
	double t = 0.5*dx.squaredNorm();
	double fp, fpp;
	switch (sepEType)
	{
		case SeparationEnergy::LOG:
			break;
		case SeparationEnergy::QUADRATIC:
			break;
		case SeparationEnergy::FLAT_LOG:
		{
			flat_log_single_hessian(xi, xj, h);
			//fp = delta / ((t + delta) * (2*t + delta));
			//fpp = -fp * fp * (3 + ((4 * t) / delta));
			//fpp = -(delta * (3 * delta + 4 * t) / ((delta + t)*(delta + t)*(delta + 2 * t)*(delta + 2 * t)));
			break;
		}
		case SeparationEnergy::QUOTIENT:
		{
			fp = delta / ((t + delta)*(t + delta));
 			fpp = -2 * fp / (t + delta);
			Mat4 Esep4;
			Esep4 << 1, 0, -1, 0,
				     0, 1, 0, -1,
				    -1, 0, 1,  0,
				     0, -1, 0, 1;
			h = fpp*dxx*dxx.transpose() + fp*Esep4;
			h = fp*Esep4;
			break;
		}
		case SeparationEnergy::QUOTIENT_NEW:
		{
			fp = delta / ((t + delta)*(t + delta));
			Mat4 Esep4;
			Esep4 << 1, 0, -1, 0,
				0, 1, 0, -1,
				-1, 0, 1, 0,
				0, -1, 0, 1;
			h = fp*Esep4;
			break;
		}
		default:
			break;
	}
	if(sepEType!=SeparationEnergy::QUOTIENT_NEW)
 		make_spd(h);
}

void Separation::flat_log_single_hessian(const Vec2& xi, const Vec2& xj, Mat4& h)
{
	double xi1 = xi(0), xi2 = xi(1);
	double xj1 = xj(0), xj2 = xj(1);
	double t4 = xi1 - xj1;
	double t2 = abs(t4);
	double t6 = xi2 - xj2;
	double t3 = abs(t6);
	double t5 = t2*t2;
	double t7 = t3*t3;
	double t8 = delta + t5 + t7;
	double t9 = 1.0 / t8;
	double t10 = sign(t4);
	double t11 = t5 + t7;
	double t16 = 1.0 / (t8*t8);
	double t22 = t2*t9*t10*2.0;
	double t23 = t2*t10*t11*t16*2.0;
	double t12 = t22 - t23;
	double t13 = t9*t11;
	double t14 = t13 + 1.0;
	double t15 = t10*t10;
	double t17 = dirac(t4);
	double t18 = 1.0 / t14;
	double t19 = sign(t6);
	double t20 = 1.0 / (t8*t8*t8);
	double t21 = 1.0 / (t14*t14);
	double t24 = t12*t12;
	double t25 = t9*t15*2.0;
	double t26 = t2*t9*t17*4.0;
	double t27 = t5*t11*t15*t20*8.0;
	double t48 = t11*t15*t16*2.0;
	double t49 = t5*t15*t16*8.0;
	double t50 = t2*t11*t16*t17*4.0;
	double t28 = t25 + t26 + t27 - t48 - t49 - t50;
	double t29 = Lsep*t18*t28;
	double t30 = t2*t3*t10*t16*t19*8.0;
	double t34 = t2*t3*t10*t11*t19*t20*8.0;
	double t31 = t30 - t34;
	double t32 = t3*t9*t19*2.0;
	double t36 = t3*t11*t16*t19*2.0;
	double t33 = t32 - t36;
	double t35 = Lsep*t18*t31;
	double t37 = Lsep*t12*t21*t33;
	double t38 = t19*t19;
	double t39 = dirac(t6);
	double t40 = t35 + t37;
	double t41 = t33*t33;
	double t42 = t9*t38*2.0;
	double t43 = t3*t9*t39*4.0;
	double t44 = t7*t11*t20*t38*8.0;
	double t54 = t11*t16*t38*2.0;
	double t55 = t7*t16*t38*8.0;
	double t56 = t3*t11*t16*t39*4.0;
	double t45 = t42 + t43 + t44 - t54 - t55 - t56;
	double t46 = Lsep*t18*t45;
	double t47 = Lsep*t21*t24;
	double t51 = -t29 + t47;
	double t52 = -t35 - t37;
	double t53 = Lsep*t21*t41;
	double t57 = -t46 + t53;
	h << t29 - Lsep*t21*t24, t52, t51, t40, -Lsep*t18*t31 - Lsep*t12*t21*t33, t46 - Lsep*t21*t41, t40, t57, t51, t40, t29 - t47, t52, t40, t57, t52, t46 - t53;
}

inline int Separation::sign(double val)
{
	return (0.0f < val) - (val < 0.0f);
}

inline double Separation::dirac(double val)
{
	return val == 0 ? INF : 0;
}

void Separation::make_spd(Mat4& h)
{
	Eigen::SelfAdjointEigenSolver<Mat4> es(h, Eigen::EigenvaluesOnly);
	Vec4 D = es.eigenvalues();
	double min_ev = D.minCoeff();
	if (min_ev < 0)
		h -= Mat4::Identity()*(min_ev - 1e-6);
}

void Separation::add_to_global_hessian(const Mat4& sh, int idx_xi, int idx_xj, int n, list<Tripletd>& htriplets)
{
	// do column by column of single-face hessian sh
	htriplets.push_back(Tripletd(idx_xi,			idx_xi,			sh(0, 0)));
	htriplets.push_back(Tripletd(idx_xi + n,	idx_xi,			sh(1, 0)));
	htriplets.push_back(Tripletd(idx_xj,			idx_xi,			sh(2, 0)));
	htriplets.push_back(Tripletd(idx_xj + n,	idx_xi,			sh(3, 0)));

	htriplets.push_back(Tripletd(idx_xi,			idx_xi + n, sh(0, 1)));
	htriplets.push_back(Tripletd(idx_xi + n,	idx_xi + n, sh(1, 1)));
	htriplets.push_back(Tripletd(idx_xj,			idx_xi + n, sh(2, 1)));
	htriplets.push_back(Tripletd(idx_xj + n,	idx_xi + n, sh(3, 1)));

	htriplets.push_back(Tripletd(idx_xi,			idx_xj,			sh(0, 2)));
	htriplets.push_back(Tripletd(idx_xi + n,	idx_xj,			sh(1, 2)));
	htriplets.push_back(Tripletd(idx_xj,			idx_xj,			sh(2, 2)));
	htriplets.push_back(Tripletd(idx_xj + n,	idx_xj,			sh(3, 2)));

	htriplets.push_back(Tripletd(idx_xi,			idx_xj + n, sh(0, 3)));
	htriplets.push_back(Tripletd(idx_xi + n,	idx_xj + n, sh(1, 3)));
	htriplets.push_back(Tripletd(idx_xj,			idx_xj + n, sh(2, 3)));
	htriplets.push_back(Tripletd(idx_xj + n,	idx_xj + n, sh(3, 3)));
}

void Separation::update_alphas(const Mat& weights, double max_possible)
{
	// factor of strengthening the separation with alpha * ||xi-xj||^2
	connect_alphas = Esep.cwiseAbs()*weights.col(0);

	// the disconnect_alphas do appear as factors for the whole
	// energy, same as edge_lengths_per_pair values
	// so it makes sense to scale them to [min, 1].
	// whereas 1 is the default, i.e. uncolored/uneffected
	// case. for every colored corner, the energy is reduced
	// and thus its min =< disconnect_alpha < 1.
	// the input comes in the range of [0, max], where max <= max_possible
	// we want to map 0 -> 1 and max_possible -> 0
	// since a factor of 0 would free the triangle we only scale
	// to [s, 1] instead, by scaling max_possible up a bit
	disconnect_alphas = 1. - (Esep.cwiseAbs()*weights.col(1) / (2. * 1.1 * max_possible)).array();
}