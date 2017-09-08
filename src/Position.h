#pragma once

#ifndef POSITION_H
#define POSITION_H

#include "EigenTypes.h"
#include "autodiff.h"

#include <map>
using namespace std;

typedef Vec Gradient;
typedef Mat Hessian;
typedef DScalar2<double, Gradient, Hessian> DScalar;

class Position
{
public:
	enum class eval_mode { F = 0, FG, FGH };

	Position();

	void init(const MatX3i& Fs_, int n);
	void value(const MatX2& X, double& f);
	void gradient(const MatX2& X, Vec& g);
	void hessian(const MatX2& X, SpMat& h);

	void evaluate_fgh(const MatX2& X, double& f, Vec& g, SpMat& h, eval_mode mode);

	int active_triangle = -1;
	Mat32 target_pos;
	map<int, Mat32> active_triangles;
	map<int, Mat32> fixed_triangles;

	double max_el = 0.0;

	vector<int> II, JJ;
	vector<double> SS;

private:
	void evaluate_single(const Mat32& current_pos, const Mat32& target_pos, double&f, Vec& g, Mat& h, eval_mode mode);

	void add_to_gradient(const Vec& gs, Vec& g, const RVec3i& verts);
	void add_to_hessian(const Mat& hs, vector<Tripletd>& htriplets, const RVec3i& verts);
	
	void prepare_hessian(int n);
	void put_into_SS(const Mat6& hs, int face);
	void process_triangle_map(const map<int, Mat32>& triangles, const MatX2& X, double& f, Vec& g, eval_mode mode);

	MatX3i Fs;
	RVec3i verts;

	Mat32 current_pos;
	RVec2 current_bary, target_bary, diff_bary;
	vector<RVec2> fixed_diff_baries;
	vector<RVec3i> vec_fixed_verts;

	int n;
	Mat6 const_single_hessian;
};

#endif