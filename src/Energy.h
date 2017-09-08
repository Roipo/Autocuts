#pragma once

#include "BBox.h"
#include "Position.h"
#include "Separation.h"
#include "EigenTypes.h"
#include "EnergySymDir.h"

#include <memory>

using namespace std;

class Energy
{
public:
	Energy();

	void init(unsigned int nf, const MatX2& Vs, const MatX3i& Fs, const MatX3& V, const MatX3i& F);
	void evaluate_f(const Vec& x, double& f);
	void evaluate_fgh(const Vec& x, double& f, Vec& g, SpMat& h);

	// helper functions
	inline void map_to_X(const Vec& x);

	unique_ptr<Separation> separation;
	unique_ptr<Position> position;
	unique_ptr<DistortionSymDir> symDirichlet;
	unique_ptr<BBox> bbox;

	// Internally work with two-column matrix instead
	// of a vector, which is used in the solver
	MatX2 X;

	double lambda = 0.0;
	double pos_weight = 100.0;
	double bbox_weight = 100.0;

	// Information on certain energy measures
	double max_sep = 0.0, max_dist = 0.0, max_pos = 0.0, max_bbox = 0.;
	double grad_norm_sep = 0.0, grad_norm_dist = 0.0, grad_norm_pos = 0.0, grad_norm_bbox = 0.;
};