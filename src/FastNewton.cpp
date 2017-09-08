#include "FastNewton.h"
#include <chrono>
#include <igl\flip_avoiding_line_search.h>



FastNewton::FastNewton()
{
}

int FastNewton::step()
{
	
	eval_fgh(m_x, f, g, h);
	//solver.compute(h);
	//p = solver.solve(-g);
	//return solver.info();
	extract_ss_from_matrix(h.triangularView<Eigen::Upper>(), SS);
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
	return 0;
}

bool FastNewton::test_progress()
{
	return true;//g.norm() > 1e-6 && diff_norm > 1e-6;
}

void FastNewton::internal_init()
{
	pardiso = make_unique<PardisoSolver<vector<int>, vector<double>>>();
	// 2 <==> real and symmetric positive definite
	// (-2 <==> real and symmetric indefinite)
	pardiso->set_type(2, true);

	// find pattern for initialization
	eval_fgh(m_x, f, g, h);
	extract_ijss_from_matrix(h.triangularView<Eigen::Upper>(), II, JJ, SS);
	pardiso->set_pattern(II, JJ, SS);
	pardiso->analyze_pattern();
}
void FastNewton::generate_triplets()
{

}
void FastNewton::internal_update_external_mesh()
{
	diff_norm = (ext_x - m_x).norm();
	ext_x = m_x;
}

void FastNewton::linesearch()
{
	Mat m_x2 = Eigen::Map<MatX2>(m_x.data(), m_x.rows() / 2, 2);
	Mat p2 = Eigen::Map<const MatX2>(p.data(), p.rows() / 2, 2);
	Mat m_plus_p = m_x2 + p2;
	double alpha = igl::flip_avoiding_line_search(Fs, m_x2, m_plus_p, bind(&FastNewton::eval_ls, this, placeholders::_1));
	m_x = Eigen::Map<Vec>(m_x2.data(), m_x2.rows() * m_x2.cols());
}

double FastNewton::eval_ls(Mat& x)
{
	double f;
	Vec g;
	SpMat h;
	Vec vec_x = Eigen::Map<Vec>(x.data(), x.rows()  * x.cols(), 1);
	eval_f(vec_x, f);
	return f;
}