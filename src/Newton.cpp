#include "Newton.h"

#include <chrono>
#include <igl/flip_avoiding_line_search.h>

Newton::Newton() {}

int Newton::step()
{
	eval_fgh(m_x, f, g, h);

 	SSd = energy->symDirichlet->SS;
 	mult(SSd, 1. - energy->lambda);

	SSs = energy->separation->SS;
	mult(SSs, energy->lambda);

	SSp = energy->position->SS;
	mult(SSp, energy->pos_weight);

	SSb = energy->bbox->SS;
	mult(SSb, energy->bbox_weight);

	SS.clear();
	SS.insert(SS.end(), SSd.begin(), SSd.end());
	SS.insert(SS.end(), SSs.begin(), SSs.end());
	SS.insert(SS.end(), SSp.begin(), SSp.end());
	SS.insert(SS.end(), SSb.begin(), SSb.end());

#ifdef USE_PARDISO
	pardiso->update_a(SS);
	try
	{
		pardiso->factorize();
	}
	catch (runtime_error& err)
	{
		cout << err.what();
		return -1;
	}
	Vec rhs = -g;
	pardiso->solve(rhs, p);
#else
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(II.size());
	int rows = *std::max_element(II.begin(), II.end())+1;
	int cols = *std::max_element(JJ.begin(), JJ.end())+1;
	assert(rows == cols && "Rows == Cols at Newton internal init");
	for(int i=0; i<II.size(); i++)
		tripletList.push_back(T(II[i],JJ[i],SS[i]));
	SpMat mat(rows, cols);
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	solver.factorize(mat);
	Vec rhs = -g;
	p = solver.solve(rhs);
#endif
	return 0;
}

bool Newton::test_progress()
{
	return true;//g.norm() > 1e-6 && diff_norm > 1e-6;
}

void Newton::internal_init()
{
#ifdef USE_PARDISO
	bool needs_init = pardiso == nullptr;

	if (needs_init)
	{
		pardiso = make_unique<PardisoSolver<vector<int>, vector<double>>>();
		pardiso->set_type(2, true);
	}
#endif

	eval_fgh(m_x, f, g, h);

	IId = energy->symDirichlet->II;
	JJd = energy->symDirichlet->JJ;
	SSd = energy->symDirichlet->SS;

	IIs = energy->separation->II;
	JJs = energy->separation->JJ;
	SSs = energy->separation->SS;

	IIp = energy->position->II;
	JJp = energy->position->JJ;
	SSp = energy->position->SS;

	IIb = energy->bbox->II;
	JJb = energy->bbox->JJ;
	SSb = energy->bbox->SS;

	if (needs_init)
	{ 
		// find pattern for initialization
		II.insert(II.end(), IId.begin(), IId.end());
		II.insert(II.end(), IIs.begin(), IIs.end());
		II.insert(II.end(), IIp.begin(), IIp.end());
		II.insert(II.end(), IIb.begin(), IIb.end());

		JJ.insert(JJ.end(), JJd.begin(), JJd.end());
		JJ.insert(JJ.end(), JJs.begin(), JJs.end());
		JJ.insert(JJ.end(), JJp.begin(), JJp.end());
		JJ.insert(JJ.end(), JJb.begin(), JJb.end());

		SS.insert(SS.end(), SSd.begin(), SSd.end());
		SS.insert(SS.end(), SSs.begin(), SSs.end());
		SS.insert(SS.end(), SSp.begin(), SSp.end());
		SS.insert(SS.end(), SSb.begin(), SSb.end());
#ifdef USE_PARDISO
		pardiso->set_pattern(II, JJ, SS);
		pardiso->analyze_pattern();
#else
		typedef Eigen::Triplet<double> T;
		std::vector<T> tripletList;
		tripletList.reserve(II.size());
		int rows = *std::max_element(II.begin(), II.end()) + 1;
		int cols = *std::max_element(JJ.begin(), JJ.end()) + 1;
		assert(rows == cols && "Rows == Cols at Newton internal init");
		for(int i=0; i<II.size(); i++)
			tripletList.push_back(T(II[i],JJ[i],SS[i]));
		SpMat mat(rows, cols);
		mat.setFromTriplets(tripletList.begin(), tripletList.end());
		solver.analyzePattern(mat);
		needs_init = false;
#endif
	}
}

void Newton::internal_update_external_mesh()
{
	diff_norm = (ext_x - m_x).norm();
	ext_x = m_x;
}

void Newton::linesearch()
{
	Mat m_x2 = Eigen::Map<MatX2>(m_x.data(), m_x.rows() / 2, 2);
	Mat p2 = Eigen::Map<const MatX2>(p.data(), p.rows() / 2, 2);
	Mat m_plus_p = m_x2 + p2;
	double alpha = igl::flip_avoiding_line_search(Fs, m_x2, m_plus_p, bind(&Newton::eval_ls, this, placeholders::_1));
	m_x = Eigen::Map<Vec>(m_x2.data(), m_x2.rows() * m_x2.cols());
}

double Newton::eval_ls(Mat& x)
{
	double f;
	Vec g;
	SpMat h;
	Vec vec_x = Eigen::Map<Vec>(x.data(), x.rows()  * x.cols(), 1);
	eval_f(vec_x, f);
	return f;
}

void Newton::mult(vector<double>& v, double s)
{
	for (int i = 0; i < v.size(); ++i)
		v[i] *= s;
}