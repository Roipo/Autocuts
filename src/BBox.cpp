#include <memory>
#include "BBox.h"

void BBox::init(int n)
{
	num_verts = n;
	box_assignments = vector<int>(num_verts, -1);
	prepare_hessian(num_verts);
}

void BBox::add_box(const MatX2& X, double xmin, double xmax, double ymin, double ymax)
{
	boxes.push_back(std::make_unique<Box>(xmin, xmax, ymin, ymax));
	int id = boxes.size() - 1;

	// assign all contained vertices of the new box to it
	for (int i = 0; i < X.rows(); ++i)
	{
		RVec vert = X.row(i);

		double xlim, ylim;
		double fx = eval_f(X(i, 0), xmin, xmax, xlim);
		double fy = eval_f(X(i, 1), ymin, ymax, ylim);

		if (fx == 0. && fy == 0.)
			box_assignments[i] = id;
	}

	box_exists = true;
}

void BBox::update_box(int id, double xmin, double xmax, double ymin, double ymax)
{
	boxes[id]->xmin = xmin;
	boxes[id]->xmax = xmax;
	boxes[id]->ymin = ymin;
	boxes[id]->ymax = ymax;
}

double BBox::eval_f(double val, double min, double max, double& lim)
{
	// in bounding box
	if (val >= min && val <= max)
	{
		lim = 0.;
		return 0.;
	}

	// below minimum
	if (val < min)
	{
		lim = min;
		return pow(min - val, 2.0);
	}

	// above maximum
	if (val > max)
	{
		lim = max;
		return pow(val - max, 2.0);
	}

	// invalid region, we should never reach this
	return INF;
}

void BBox::value(const MatX2& X, double& f)
{
	f = 0.0;

	if (!box_exists)
		return;

	missed_limit = MatX2::Zero(X.rows(), X.cols());

	int threads = omp_get_max_threads();
	Vec f_per_thread = Vec::Zero(threads);
	Vec fmax_per_thread = Vec::Zero(threads);

	#pragma omp parallel for num_threads(threads)
	for (int i = 0; i < X.rows(); ++i)
	{
		int tid = omp_get_thread_num();

		RVec vert = X.row(i);

		int assignment = box_assignments[i];
		if (assignment == -1)
			continue;

		Box* b = boxes[assignment].get();

		double xlim, ylim;
		double fx = eval_f(X(i, 0), b->xmin, b->xmax, xlim);
		double fy = eval_f(X(i, 1), b->ymin, b->ymax, ylim);

		missed_limit(i, 0) = xlim;
		missed_limit(i, 1) = ylim;

		double fi = fx + fy;
		f_per_thread(tid) +=  fi;

		if (fmax_per_thread(tid) < fi)
			fmax_per_thread(tid) = fi;
	}
	f_max = fmax_per_thread.maxCoeff();
	f = f_per_thread.sum();
}

void BBox::gradient(const MatX2& X, Vec& g)
{
	if (!box_exists)
	{
		g = Vec::Zero(X.rows() * 2);
		return;
	}

	MatX2 ge = MatX2::Zero(X.rows(), X.cols());
	int threads = omp_get_max_threads();
	#pragma omp parallel for num_threads(threads)
	for (int i = 0; i < X.rows(); ++i)
	{
		RVec2 out = missed_limit.row(i);
		RVec vert = X.row(i);

		int assignment = box_assignments[i];
		if (assignment == -1)
			continue;

		Box* b = boxes[assignment].get();

		if (out(0))
			ge(i, 0) = 2 * (vert(0) - out(0));
		if (out(1))
			ge(i, 1) = 2 * (vert(1) - out(1));
	}
	g = Eigen::Map<Vec>(ge.data(), ge.rows() * 2, 1);
}

void BBox::hessian(const MatX2& X)
{
	SS = vector<double>(II.size(), 0.);

	if (!box_exists)
		return;

	int n = X.rows();
	int threads = omp_get_max_threads();
#pragma omp parallel for num_threads(threads)
	for (int i = 0; i < X.rows(); ++i)
	{
		RVec2 out = missed_limit.row(i);

		if (out(0))
			SS[i] = 2.;
		if (out(1))
			SS[i + n] = 2.;
	}
}

void BBox::prepare_hessian(int n)
{
	II.clear();
	JJ.clear();
	auto PushPair = [&](int i, int j) { II.push_back(i); JJ.push_back(j); };
	for (int i = 0; i < 2 * n; ++i)
	{
		// Since each vertex is treated solely on its own
		// the hessian is simply a diagonal matrix with
		// either a 0 if the vertex/coordinate is in the box
		// or a 2, because d/dx^2 x^2 = d/dx 2x = 2
		PushPair(i, i);
	}
	SS = vector<double>(II.size(), 0.);
}