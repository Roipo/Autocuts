#include "Position.h"

#include <list>
#include <igl/slice.h>
#include <igl/slice_into.h>

DECLARE_DIFFSCALAR_BASE();

Position::Position()
{
	const double v = 2.0 / 9.0;
	const_single_hessian << v, v, v, 0, 0, 0,
		v, v, v, 0, 0, 0,
		v, v, v, 0, 0, 0,
		0, 0, 0, v, v, v,
		0, 0, 0, v, v, v,
		0, 0, 0, v, v, v;
	const_single_hessian += Mat6::Identity()*1e-4;
}

void Position::init(const MatX3i& Fs_, int n)
{
	Fs = Fs_;
	prepare_hessian(n);
}

void Position::value(const MatX2& X, double& f)
{
	f = 0.0;
	if (active_triangle == -1 && fixed_triangles.size() == 0)
		return;

	if (active_triangle > -1)
	{
		// treat active triangle
		verts = Fs.row(active_triangle);
		igl::slice(X, verts, 1, current_pos);

		target_bary = target_pos.colwise().mean();
		current_bary = current_pos.colwise().mean();
		diff_bary = current_bary - target_bary;

		f += (diff_bary).squaredNorm();
	}

	// add fixed triangles
	fixed_diff_baries.clear();
	vec_fixed_verts.clear();
	RVec2 fixed_diff_bary;
	RVec3i fixed_verts;
	for (pair<int, Mat32> p : fixed_triangles)
	{
		fixed_verts = Fs.row(p.first);
		igl::slice(X, fixed_verts, 1, current_pos);
		vec_fixed_verts.push_back(fixed_verts);

		target_bary = p.second.colwise().mean();
		current_bary = current_pos.colwise().mean();
		fixed_diff_bary = current_bary - target_bary;

		f += (fixed_diff_bary).squaredNorm();
		fixed_diff_baries.push_back(fixed_diff_bary);
	}
}

void Position::gradient(const MatX2& X, Vec& g)
{
	g = Vec::Zero(2 * X.rows());
	if (active_triangle == -1 && fixed_triangles.size() == 0)
		return;

	MatX2 gp = MatX2::Zero(X.rows(), 2);

	RVec2 grad;
	if (active_triangle > -1)
	{
		// treat active triangle
		grad = 2.0 / 3.0 * diff_bary;
		gp.row(verts(0)) = grad;
		gp.row(verts(1)) = grad;
		gp.row(verts(2)) = grad;
	}

	// treat fixed triangles
	RVec2 fixed_diff_bc;
	RVec3i fixed_verts;
	int i = 0;
	for (pair<int, Mat32> p : fixed_triangles)
	{
		fixed_diff_bc = fixed_diff_baries.at(i);
		fixed_verts = vec_fixed_verts.at(i);

		grad = 2.0 / 3.0 * fixed_diff_bc;
		gp.row(fixed_verts(0)) += grad;
		gp.row(fixed_verts(1)) += grad;
		gp.row(fixed_verts(2)) += grad;

		++i;
	}

	g = Eigen::Map<Vec>(gp.data(), gp.rows() * 2, 1);
}

void Position::hessian(const MatX2& X, SpMat& h)
{
	vector<Tripletd> htriplets;

	// all hessian values are 2/9, the ony task is to distribute them correctly
	
	if (active_triangle > -1)
	{ // treat active triangle
		int add = 0;
		int addl = 0;
		for (int k = 0; k < 2; ++k)
		{
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					htriplets.push_back(Tripletd(verts(i) + add, verts(j) + add, const_single_hessian(i + addl, j + addl)));
				}
			}
			add += n;
			addl += 3;
		}
	}

	// treat fixed triangles
	RVec2 fixed_diff_bc;
	RVec3i fixed_verts;
	int i = 0;
	for (pair<int, Mat32> p : fixed_triangles)
	{
		fixed_verts = vec_fixed_verts.at(i);

		int add = 0;
		int addl = 0;
		for (int k = 0; k < 2; ++k)
		{
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					htriplets.push_back(Tripletd(verts(i) + add, verts(j) + add, const_single_hessian(i + addl, j + addl)));
				}
			}
			add += n;
			addl += 3;
		}

		++i;
	}
	h = SpMat(2 * X.rows(), 2 * X.rows());
	h.setFromTriplets(htriplets.begin(), htriplets.end());
}

void Position::evaluate_fgh(const MatX2& X, double& f, Vec& g, SpMat& h, eval_mode mode)
{
	n = X.rows();
	f = 0.0;
	g = Vec::Zero(2 *n);
	max_el = 0.0;

	SS = vector<double>(II.size(), 0.);
	if (active_triangle == -1 && fixed_triangles.size() == 0)
		return;

	double fs;
	Vec gs;
	Mat hs;

	DiffScalarBase::setVariableCount(6);

	if (active_triangle > -1)
	{
		verts = Fs.row(active_triangle);
		igl::slice(X, verts, 1, current_pos);

		evaluate_single(current_pos, target_pos, fs, gs, hs, mode);

		f += fs;
		if (fs > max_el)
			max_el = fs;

		if (mode >= eval_mode::FG)
				add_to_gradient(gs, g, verts);
		if (mode == eval_mode::FGH)
		{
			Mat6 hsl = hs;
			hsl += Mat6::Identity()*1e-4;
			put_into_SS(hsl, active_triangle);
		}
	}
	process_triangle_map(fixed_triangles, X, f, g, mode);
}

void Position::evaluate_single(const Mat32& current_pos, const Mat32& target_pos, double&f, Vec& g, Mat& h, eval_mode mode)
{
	DScalar p11(0, current_pos(0, 0)), p21(1, current_pos(1, 0)), p31(2, current_pos(2, 0));
	DScalar p12(3, current_pos(0, 1)), p22(4, current_pos(1, 1)), p32(5, current_pos(2, 1));

	double t11 = target_pos(0, 0), t21 = target_pos(1, 0), t31 = target_pos(2, 0);
	double t12 = target_pos(0, 1), t22 = target_pos(1, 1), t32 = target_pos(2, 1);

	DScalar current_bary_x = (p11 + p21 + p31) / 3, current_bary_y = (p12 + p22 + p32) / 3;
	double target_bary_x = (t11 + t21 + t31) / 3, target_bary_y = (t12 + t22 + t32) / 3;

	DScalar diff_bary_x = current_bary_x - target_bary_x;
	DScalar diff_bary_y = current_bary_y - target_bary_y;

	DScalar fx = diff_bary_x*diff_bary_x + diff_bary_y*diff_bary_y;

	f = fx.getValue();
	if (mode >= eval_mode::FG)
		g = fx.getGradient();
	if (mode == eval_mode::FGH)
		h = fx.getHessian();
}

void Position::add_to_gradient(const Vec& gs, Vec& g, const RVec3i& verts)
{
	g(verts(0)) += gs(0);
	g(verts(1)) += gs(1);
	g(verts(2)) += gs(2);
	g(verts(0) + n) += gs(3);
	g(verts(1) + n) += gs(4);
	g(verts(2) + n) += gs(5);
}

void Position::add_to_hessian(const Mat& hs, vector<Tripletd>& htriplets, const RVec3i& verts)
{
	int add = 0;
	int addl = 0;
	for (int k = 0; k < 2; ++k)
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				htriplets.push_back(Tripletd(verts(i) + add, verts(j) + add, hs(i + addl, j + addl)));
			}
		}
		add += n;
		addl += 3;
	}
}

void Position::prepare_hessian(int n)
{
	II.clear();
	JJ.clear();
	auto PushPair = [&](int i, int j) { II.push_back(i); JJ.push_back(j); };
	for (int i = 0; i < Fs.rows(); ++i)
	{
		// the single face hessians are 6x6 block diagonal matrices with
		// two 3x3 blocks (filled with constants 0.222222 == 2/9)
		// since we only need the upper diagonal of these 6x6, respective
		// two 3x3 blocks only 12 values are needed. they are accessed
		// and also put into the big hessian in column order

		// base indices, x and y indicate second derivate w.r.t x or y
		// mixed second derivates d/dxy = d/dyx = 0
		int xbi = 3 * i; // x base index
		int ybi = 3 * i + n; // y base index

		// first column
		PushPair(xbi,			xbi);

		// second column
		PushPair(xbi,			xbi + 1);
		PushPair(xbi + 1, xbi + 1);

		// third column
		PushPair(xbi,			xbi + 2);
		PushPair(xbi + 1, xbi + 2);
		PushPair(xbi + 2, xbi + 2);

		// fourth column
		PushPair(ybi,			ybi);

		// fifth column
		PushPair(ybi,			ybi + 1);
		PushPair(ybi + 1, ybi + 1);

		// sixth column
		PushPair(ybi,			ybi + 2);
		PushPair(ybi + 1, ybi + 2);
		PushPair(ybi + 2, ybi + 2);
	}
	SS = vector<double>(II.size(), 0.);
}

void Position::put_into_SS(const Mat6& hs, int face)
{
	int index = 12 * face;
	// upper-left (x) block
	for (int a = 0; a < 3; ++a)
	{
		for (int b = 0; b <= a; ++b)
		{
			SS[index++] = hs(b, a);
		}
	}
	// lower-right (y) block
	for (int a = 3; a < 6; ++a)
	{
		for (int b = 3; b <= a; ++b)
		{
			SS[index++] = hs(b, a);
		}
	}
}

void Position::process_triangle_map(const map<int, Mat32>& triangles, const MatX2& X, double& f, Vec& g, eval_mode mode)
{
	double fs;
	Vec gs;
	Mat hs;
	RVec3i verts;
	Mat32 tar;
	for (pair<int, Mat32> p : triangles)
	{
		verts = Fs.row(p.first);
		igl::slice(X, verts, 1, current_pos);

		tar = p.second;
		evaluate_single(current_pos, tar, fs, gs, hs, mode);

		f += fs;
		if (fs > max_el)
			max_el = fs;

		if (mode >= eval_mode::FG)
			add_to_gradient(gs, g, verts);
		if (mode == eval_mode::FGH)
		{
			Mat6 hsl = hs;
			hsl += Mat6::Identity()*1e-4;
			//add_to_hessian(hsl, htriplets, fixed_verts);
			put_into_SS(hsl, p.first);
		}
	}
}