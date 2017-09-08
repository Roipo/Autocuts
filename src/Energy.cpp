#include "Energy.h"

Energy::Energy()
	:
	separation(make_unique<Separation>()),
	position(make_unique<Position>()),
	symDirichlet(make_unique<DistortionSymDir>()),
	bbox(make_unique<BBox>())
{
}

void Energy::init(unsigned int nf, const MatX2& Vs, const MatX3i& Fs, const MatX3& V, const MatX3i& F)
{
	separation->init(Vs.rows());
	position->init(Fs, Vs.rows());
	symDirichlet->init(V, F, Vs, Fs);
	bbox->init(Vs.rows());
}

void Energy::evaluate_f(const Vec& x, double& f)
{
	map_to_X(x);

	double fs;
	separation->value(X, fs);

	double fd;
	Vec gd;
	SpMat hd;
	symDirichlet->value(X, fd);

	double fp;
	Vec gp;
	SpMat hp;
	position->evaluate_fgh(X, fp, gp, hp, Position::eval_mode::F);

	double fb;
	bbox->value(X, fb);

	f = (1.0 - lambda)*fd + lambda*fs + pos_weight*fp + bbox_weight*fb;
}

void Energy::evaluate_fgh(const Vec& x, double& f, Vec& g, SpMat& h)
{
	map_to_X(x);

	double fs, fd;
	Vec gs, gd;
	SpMat hd;
	
	// Separation
	separation->value(X, fs);
	separation->gradient(X, gs);
	separation->hessian(X);

	max_sep = separation->f_per_pair.maxCoeff();
	grad_norm_sep = gs.norm();

	// Distortion
	symDirichlet->value(X, fd);
	symDirichlet->gradient(X, gd);
	symDirichlet->hessian(X);

	max_dist = symDirichlet->Efi.maxCoeff();
	grad_norm_dist = gd.norm();

	// Position
	double fp;
	Vec gp;
	SpMat hp;
	position->evaluate_fgh(X, fp, gp, hp, Position::eval_mode::FGH);
	max_pos = position->max_el;
	grad_norm_pos = gp.norm();

	// bbox
	double fb;
	Vec gb;
	bbox->value(X, fb);
	bbox->gradient(X, gb);
	bbox->hessian(X);
	max_bbox = bbox->f_max;
	grad_norm_bbox = gb.norm();

	f = (1.0 - lambda)*fd + lambda*fs + pos_weight*fp + bbox_weight*fb;
	g = (1.0 - lambda)*gd + lambda*gs + pos_weight*gp + bbox_weight*gb;
}

inline void Energy::map_to_X(const Vec& x)
{
	X = Eigen::Map<const MatX2>(x.data(), x.rows() / 2, 2);
}