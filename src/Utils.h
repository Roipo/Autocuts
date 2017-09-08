#pragma once

#ifndef UTILS_H
#define UTILS_H

#include "EigenTypes.h"

#include <vector>
#include <igl\sum.h>
#include <igl\diag.h>
#include <igl\arap.h>
#include <igl\lscm.h>
#include <igl\grad.h>
#include <igl\slice.h>
#include <igl\unique.h>
#include <igl\sparse.h>
#include <igl\harmonic.h>
#include <igl\sortrows.h>
#include <igl\slice_into.h>
#include <igl\local_basis.h>
#include <igl\boundary_loop.h>
#include <igl\per_face_normals.h>
#include <igl\adjacency_matrix.h>
#include <igl\map_vertices_to_circle.h>
#include <igl\project_isometrically_to_plane.h>

// for EXCEPTION_POINTERS
#include <Windows.h>

using namespace std;

class Utils
{
public:
	enum class Init : int { RANDOM, ISOMETRIC, LOADED, LOAD_FROM_FILE };

	static void get_soup_and_helper_matrices(const MatX3& V_in, const MatX3i& F, SpMat& EVvar1, SpMat& EVvar2, SpMat& V2V)
	{
		MatX3 V, Vs3d_similar;
		MatX3i Fs;
		init_mesh(V_in, F, V, E);
		generate_soup_3d_similar(V, F, Vs3d_similar, Fs);
		compute_edges(Fs, Es);
		nvs = Vs3d_similar.rows(); nfs = Fs.rows(); nes = Es.rows();
		compute_V2V(F, V2V);
		Mat E2Edt;
		compute_E2E(V2V, E2Edt);
		compute_EVvars(Fs, E2Edt, EVvar1, EVvar2);
	}

	static void Utils::compute_EVvars(const MatX3i& Fs, const Mat& E2Edt, SpMat& EVvar1, SpMat& EVvar2)
	{
		Mati El(Fs.rows(), 6);
		El << Fs.col(0), Fs.col(1), Fs.col(2), Fs.col(0), Fs.col(1), Fs.col(2);
		Mati Elt = El.transpose();
		MatX2i Elf = Eigen::Map<Mat2Xi>(Elt.data(), 2, Elt.cols() * Elt.rows() / 2.0).transpose();

		SpMati K_V2E;
		Veci lin_twin(6 * nfs), E_twin(6 * nfs), S_twin(6 * nfs);

		Veci lin = Eigen::VectorXi::LinSpaced(3 * nfs, 0, 3 * nfs - 1);
		lin_twin << lin, lin;
		E_twin << Elf.col(0), Elf.col(1);
		S_twin << Eigen::VectorXi::Ones(3 * nfs), -Eigen::VectorXi::Ones(3 * nfs);

		igl::sparse(lin_twin, E_twin, S_twin, 3 * nfs, nvs, K_V2E);

		int s;
		vector<int> safe_cols;
		for (int i = 0; i < E2Edt.cols(); ++i) {
			s = 0;
			for (int j = 0; j < E2Edt.rows(); ++j) {
				s += E2Edt(j, i);
			}
			if (s == 0) { safe_cols.push_back(i); }
		}

		Mat E2Ec = Mat(E2Edt.rows(), safe_cols.size());
		for (int i = 0; i < safe_cols.size(); ++i) {
			E2Ec.col(i) = E2Edt.col(safe_cols[i]);
		}
		Mat E2Ect = E2Ec.transpose();

		SpMat E2Ecs = E2Ect.sparseView();

		SpMat K_V2Er = K_V2E.cast<double>() * E2Ect.colwise().sum().asDiagonal();
		SpMat K_V2Es = K_V2Er;
		SpMat K_V2Ed = -K_V2Er;

		K_V2Es = K_V2Es.unaryExpr([](double v) {
			return  v < 0 ? 0 : v;
		});
		K_V2Ed = K_V2Ed.unaryExpr([](double v) {
			return v < 0 ? 0 : v;
		});

		K_V2Es.prune(0.01);
		K_V2Ed.prune(0.01);

		EVvar1 = E2Ecs * K_V2Es.transpose();
		EVvar2 = E2Ecs * K_V2Ed.transpose();
	}

	static void Utils::compute_V2V(const MatX3i& F, SpMat& V2V)
	{
		Veci lin = Eigen::VectorXi::LinSpaced(3 * nf, 0, 3 * nf - 1);
		Mat3Xi Ft = F.transpose();
		Veci Fv = Eigen::Map<Veci>(Ft.data(), 3 * F.rows(), 1);
		igl::sparse(Fv, lin, Eigen::VectorXi::Ones(3 * F.rows()), V2V);
	}

	static void Utils::init_mesh(const MatX3& V_in, const MatX3i& F, MatX3& V_out, MatX2& E)
	{
		centralize_vertices(V_in, V_out);
		nv = V_out.rows();
		nf = F.rows();
		compute_edges(F, E);
		ne = E.rows();
	}

	static void Utils::generate_soup_3d_similar(const MatX3& V, const MatX3i& F, MatX3& Vs, MatX3i& Fs)
	{
		Veci lin = Eigen::VectorXi::LinSpaced(3 * nf, 0, 3 * nf - 1);
		Mat3Xi Ft = F.transpose();
		Veci Fv = Eigen::Map<Veci>(Ft.data(), 3 * nf, 1);
		igl::slice(V, Fv, 1, Vs);
		Fs = Eigen::Map<Mat3Xi>(lin.data(), 3, nf).transpose();
		nvs = Vs.rows(); nfs = Fs.rows();
	}

	static void Utils::generate_soup_2d_random(const MatX3& V, const MatX3i& F, MatX2& Vs, MatX3i& Fs)
	{
		Veci lin = Eigen::VectorXi::LinSpaced(3 * nf, 0, 3 * nf - 1);
		Fs = Eigen::Map<Mat3Xi>(lin.data(), 3, nf).transpose();

		Vs = (MatX2::Random(3 * nf, 2).array() * 2.0).array();
	}

	static void Utils::generate_soup_2d_iso(const MatX3& V, const MatX3i& F, MatX2& Vs, MatX3i& Fs)
	{
		SpMat I;
		MatX2 Vst;
		MatX3i Fst;
		igl::project_isometrically_to_plane(V, F, Vst, Fst, I);

		Vs = MatX2(Vst.rows(), 2);

		Veci lin = Eigen::VectorXi::LinSpaced(3 * nf, 0, 3 * nf - 1);
		Fs = Eigen::Map<Mat3Xi>(lin.data(), 3, nf).transpose();

		Mat32 face;
		for (int i = 0; i < Fst.rows(); ++i)
		{
			igl::slice(Vst, Fst.row(i).eval(), 1, face);
			Vs.block<3, 2>(3 * i, 0) = face;
		}
		Vs += MatX2::Random(Vs.rows(), Vs.cols())*1e-6;
	}

	static void Utils::generate_soup_2d_harmonic(const MatX3& V, const MatX3i& F, MatX2& Vs, MatX3i& Fs)
	{
		Veci bnd;
		igl::boundary_loop(F, bnd);
		Mat bnd_uv;
		igl::map_vertices_to_circle(V, bnd, bnd_uv);

		Mat initial_guess;
		igl::harmonic(V, F, bnd, bnd_uv, 1, initial_guess);

		Veci lin = Eigen::VectorXi::LinSpaced(3 * nf, 0, 3 * nf - 1);
		Fs = Eigen::Map<Mat3Xi>(lin.data(), 3, nf).transpose();

		Vs = MatX2(3 * Fs.rows(), 2);

		Mat32 face;
		for (int i = 0; i < F.rows(); ++i)
		{
			igl::slice(initial_guess, F.row(i).eval(), 1, face);
			Vs.block<3, 2>(3 * i, 0) = face;
		}
		Vs += MatX2::Random(Vs.rows(), Vs.cols())*1e-6;
	}

	static void Utils::map_vertices_to_circle_area_normalized(
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::VectorXi& bnd,
		Eigen::MatrixXd& UV) {

		Eigen::VectorXd dblArea_orig; // TODO: remove me later, waste of computations
		igl::doublearea(V, F, dblArea_orig);
		double area = dblArea_orig.sum() / 2;
		double radius = sqrt(area / (M_PI));

		// Get sorted list of boundary vertices
		std::vector<int> interior, map_ij;
		map_ij.resize(V.rows());
		interior.reserve(V.rows() - bnd.size());

		std::vector<bool> isOnBnd(V.rows(), false);
		for (int i = 0; i < bnd.size(); i++)
		{
			isOnBnd[bnd[i]] = true;
			map_ij[bnd[i]] = i;
		}

		for (int i = 0; i < (int)isOnBnd.size(); i++)
		{
			if (!isOnBnd[i])
			{
				map_ij[i] = interior.size();
				interior.push_back(i);
			}
		}

		// Map boundary to unit circle
		std::vector<double> len(bnd.size());
		len[0] = 0.;

		for (int i = 1; i < bnd.size(); i++)
		{
			len[i] = len[i - 1] + (V.row(bnd[i - 1]) - V.row(bnd[i])).norm();
		}
		double total_len = len[len.size() - 1] + (V.row(bnd[0]) - V.row(bnd[bnd.size() - 1])).norm();

		UV.resize(bnd.size(), 2);
		for (int i = 0; i < bnd.size(); i++)
		{
			double frac = len[i] * (2. * M_PI) / total_len;
			UV.row(map_ij[bnd[i]]) << radius*cos(frac), radius*sin(frac);
		}

	}

	static void Utils::seperate_3d_and_uv_mesh(const MatX3& V, const MatX3i& F, MatX2& Vs, MatX3i& Fs)
	{
		// uv mesh on the left, 3d mesh on the right
		// find right corner of uv and left of 3d mesh and seperate them by translating the uv mesh to the left
		double most_right_uv = Vs.col(0).maxCoeff();
		double most_left_3d = V.col(0).minCoeff();
		double diff = most_left_3d - most_right_uv;
		Vs.col(0).array() += 1.5*(diff + diff*0.1);
	}

	static void Utils::correct_flipped_triangles(MatX2& Vs, MatX3i& Fs)
	{
		Mat32 f_uv;
		for (int i = 0; i < Fs.rows(); ++i)
		{
			igl::slice(Vs, Fs.row(i).eval(), 1, f_uv);
			Vec2 e1uv = f_uv.row(1) - f_uv.row(0);
			Vec2 e2uv = f_uv.row(2) - f_uv.row(0);

			Vec3 uve1, uve2;
			uve1 << e1uv, 0;
			uve2 << e2uv, 0;
			Vec3 x = uve1.cross(uve2);

			if (x(2) < 0) // flipped triangle
			{
				Vs(Fs(i, 0), 1) = -Vs(Fs(i, 0), 1);
				Vs(Fs(i, 1), 1) = -Vs(Fs(i, 1), 1);
				Vs(Fs(i, 2), 1) = -Vs(Fs(i, 2), 1);
			}
		}
	}

	static void Utils::compute_E2E(const SpMat& V2V, Mat& E2Edt)
	{
		SpMati Ic_G, Ic_K;
		compute_V2E(E, Ic_G, false);
		compute_V2E(Es, Ic_K, true);
		unsigned int ic_g_c = Ic_G.cols();
		unsigned int ic_k_c = Ic_K.cols();

		// find t(i,j) == 2 and replace with #cols of Ic_G into e2ev
		SpMati t = (Ic_G.transpose() * V2V * Ic_K);
		unsigned int tr = t.rows();
		vector<int> v;
		for (int i = 0; i < t.outerSize(); ++i) {
			for (Eigen::SparseMatrix<int>::InnerIterator it(t, i); it; ++it) {
				if (it.value() == 2) {
					v.push_back((it.col() * tr + it.row() + 1) % ic_g_c);
				}
			}
		}

		// -1 needed to go back to 0-indexing
		Veci e2ev(v.size());
		for (int i = 0; i < v.size(); ++i)
			e2ev(i) = (v[i] == 0 ? ic_g_c : v[i]) - 1;

		SpMat E2E;
		Veci lin = Eigen::VectorXi::LinSpaced(ic_k_c, 0, ic_k_c - 1);
		igl::sparse(lin, e2ev, Veci::Ones(ic_k_c), ic_k_c, ic_g_c, E2E);

		// find last nnz element of each col
		unsigned int last;
		for (int i = 0; i < E2E.outerSize(); ++i) {
			for (SpMat::InnerIterator it(E2E, i); it; ++it)
			{
				last = it.row();
			}
			E2E.coeffRef(last, i) = -1;
		}
		E2Edt = E2E.toDense();
	}

	static void Utils::compute_V2E(const MatX2& E, SpMati& V2E, bool is_soup)
	{
		unsigned int ne, nv;
		if (is_soup) { ne = nes; nv = nvs; }
		else { ne = Utils::ne; nv = Utils::nv; }

		MatX2i Ei = E.cast<int>();
		Veci I(2 * ne), J(2 * ne);
		I << Ei.col(0), Ei.col(1);

		Veci lin = Eigen::VectorXi::LinSpaced(ne, 0, ne - 1);
		J << lin, lin;

		igl::sparse(I, J, Eigen::VectorXi::Ones(2 * ne), nv, ne, V2E);
	}

	static void Utils::compute_edges(const MatX3i& F, MatX2& E)
	{
		Veci I(3 * nf), J(3 * nf), IC, IA;
		MatX2i dE = MatX2i(3 * nf, 2), sE, Ei;

		// list all edge pairs (including duplicates)
		I << F.col(0), F.col(1), F.col(2);
		J << F.col(1), F.col(2), F.col(0);
		dE << I, J;

		// make E(i, 0) < E(i, 1)
		for (int i = 0; i < dE.rows(); ++i)
			if (dE(i, 0) > dE(i, 1))
				swap(dE(i, 0), dE(i, 1));

		Mati IX;
		igl::sortrows(dE, true, sE, IX);
		igl::unique_rows(sE, Ei, IA, IC);
		E = Ei.cast<double>();
	}

	static void Utils::centralize_vertices(const MatX3& V_in, MatX3& V_out)
	{
		RVec3 bary = (V_in.colwise().minCoeff() + V_in.colwise().maxCoeff()) / 2.0;
		V_out = V_in.rowwise() - bary;
		double dt = V_out.cwiseAbs().maxCoeff();
		V_out = V_out / dt;
	}

	static void Utils::computeSurfaceGradientMatrix(const MatX3& V, const MatX3i& F, SpMat& D1, SpMat& D2)
	{
		MatX3 F1, F2, F3;
		SpMat DD, Dx, Dy, Dz;

		igl::local_basis(V, F, F1, F2, F3);
		igl::grad(V, F, DD);

		Dx = DD.topLeftCorner(F.rows(), V.rows());
		Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
		Dz = DD.bottomRightCorner(F.rows(), V.rows());

		D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
		D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
	}

	static void Utils::compute_triangle_areas(const MatX3& V, const MatX3i& F, Vec& A)
	{
		Veci F1 = F.col(0), F2 = F.col(1), F3 = F.col(2);
		MatX3 VF1, VF2, VF3;
		igl::slice(V, F1, 1, VF1);
		igl::slice(V, F2, 1, VF2);
		igl::slice(V, F3, 1, VF3);
		Vec L1 = sqrt((VF2 - VF3).array().pow(2.0).rowwise().sum());
		Vec L2 = sqrt((VF1 - VF3).array().pow(2.0).rowwise().sum());
		Vec L3 = sqrt((VF1 - VF2).array().pow(2.0).rowwise().sum());
		Vec S = (L1 + L2 + L3).array() / 2.0;
		Vec SL1 = S - L1;
		Vec SL2 = S - L2;
		Vec SL3 = S - L3;
		A = sqrt(S.array() * SL1.array() * SL2.array() * SL3.array());
	}

	static void computeSurfaceGradientPerFace(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, Eigen::MatrixX3d &D1, Eigen::MatrixX3d &D2)
	{
		using namespace Eigen;
		MatrixX3d F1, F2, F3;
		igl::local_basis(V, F, F1, F2, F3);
		const int Fn = F.rows();  const int vn = V.rows();

		MatrixXd Dx(Fn, 3), Dy(Fn, 3), Dz(Fn, 3);
		MatrixXd fN; igl::per_face_normals(V, F, fN);
		VectorXd Ar; igl::doublearea(V, F, Ar);
		PermutationMatrix<3> perm;

		Vec3i Pi;
		Pi << 1, 2, 0;
		PermutationMatrix<3> P = Eigen::PermutationMatrix<3>(Pi);

		for (int i = 0; i < Fn; i++) {
			// renaming indices of vertices of triangles for convenience
			int i1 = F(i, 0);
			int i2 = F(i, 1);
			int i3 = F(i, 2);

			// #F x 3 matrices of triangle edge vectors, named after opposite vertices
			Matrix3d e;
			e.col(0) = V.row(i2) - V.row(i1);
			e.col(1) = V.row(i3) - V.row(i2);
			e.col(2) = V.row(i1) - V.row(i3);;

			Vector3d Fni = fN.row(i);
			double Ari = Ar(i);

			//grad3_3f(:,[3*i,3*i-2,3*i-1])=[0,-Fni(3), Fni(2);Fni(3),0,-Fni(1);-Fni(2),Fni(1),0]*e/(2*Ari);
			Matrix3d n_M;
			n_M << 0, -Fni(2), Fni(1), Fni(2), 0, -Fni(0), -Fni(1), Fni(0), 0;
			VectorXi R(3); R << 0, 1, 2;
			VectorXi C(3); C << 3 * i + 2, 3 * i, 3 * i + 1;
			Matrix3d res = ((1. / Ari)*(n_M*e))*P;

			Dx.row(i) = res.row(0);
			Dy.row(i) = res.row(1);
			Dz.row(i) = res.row(2);
			//		igl::slice_into(res, R, C, grad3_3f);
		}
		D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
		D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
	}

	static inline void SSVD2x2(const Eigen::Matrix2d& A, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
	{
		double e = (A(0) + A(3))*0.5;
		double f = (A(0) - A(3))*0.5;
		double g = (A(1) + A(2))*0.5;
		double h = (A(1) - A(2))*0.5;
		double q = sqrt((e*e) + (h*h));
		double r = sqrt((f*f) + (g*g));
		double a1 = atan2(g, f);
		double a2 = atan2(h, e);
		double rho = (a2 - a1)*0.5;
		double phi = (a2 + a1)*0.5;

		S(0) = q + r;
		S(1) = 0;
		S(2) = 0;
		S(3) = q - r;

		double c = cos(phi);
		double s = sin(phi);
		U(0) = c;
		U(1) = s;
		U(2) = -s;
		U(3) = c;

		c = cos(rho);
		s = sin(rho);
		V(0) = c;
		V(1) = -s;
		V(2) = s;
		V(3) = c;
	}
private:
	static unsigned int nv, nf, ne, nvs, nfs, nes;
	static MatX2 E, Es;
	static MatX3 Vs3d_similar;
};

class Timer {
	typedef std::chrono::high_resolution_clock high_resolution_clock;
	typedef std::chrono::milliseconds milliseconds;
	typedef std::chrono::microseconds microseconds;
public:
	explicit Timer(bool run = false)
	{
		if (run)
			Reset();
	}
	void Reset()
	{
		_start = high_resolution_clock::now();
	}
	microseconds Elapsed() const
	{
		return std::chrono::duration_cast<microseconds>(high_resolution_clock::now() - _start);
	}
	template <typename T, typename Traits>
	friend std::basic_ostream<T, Traits>& operator<<(std::basic_ostream<T, Traits>& out, const Timer& timer)
	{
		return out << timer.Elapsed().count();
	}
private:
	high_resolution_clock::time_point _start;
};
#endif