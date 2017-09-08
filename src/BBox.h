#pragma once

#include "EigenTypes.h"

#include <memory>

using namespace std;

#ifndef INF
#define INF numeric_limits<double>::infinity()
#endif

class BBox
{
public:
	void init(int n);
	void add_box(const MatX2& X, double xmin, double xmax, double ymin, double ymax);
	void update_box(int id, double xmin, double xmax, double ymin, double ymax);

	void value(const MatX2& X, double& f);
	void gradient(const MatX2& X, Vec& g);
	void hessian(const MatX2& X);

	void prepare_hessian(int n);

	// hessian/pardiso structures
	vector<int> II, JJ;
	vector<double> SS;

	double f_max = 0.0;

private:
	struct Box
	{
		Box(double _xmin, double _xmax, double _ymin, double _ymax)
		{
			xmin = _xmin;
			xmax = _xmax;
			ymin = _ymin;
			ymax = _ymax;
		}
		// corners
		double xmin, xmax, ymin, ymax;
	};

	double eval_f(double val, double min, double max, double& lim);

	// bookkeeping done in value function to store which vertices and
	// coordinates are out of the box.
	// if missed_limit == 0. the vertex is already in the box.
	// else (!= 0.) this is the boundary value that was missed
	// this works because x and y have the same gradient
	// 2 ('x or y value' - 2*lim)
	MatX2 missed_limit;

	// index list for each vertex to which box it belongs and thus
	// with which boundaries the computations are made
	vector<int> box_assignments;

	// a list of boxes
	vector<unique_ptr<Box>> boxes;

	int num_verts = 0;
	bool box_exists = false;
};