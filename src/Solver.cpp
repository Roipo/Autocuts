#include "Solver.h"

#include "Utils.h"

#include <iostream>
#include <igl\readOBJ.h>
#include <igl\file_dialog_open.h>
#include <igl\project_isometrically_to_plane.h>

Solver::Solver()
	:
	energy(make_shared<Energy>()),
	param_mutex(make_unique<mutex>()),
	mesh_mutex(make_unique<shared_mutex>()),
	param_cv(make_unique<condition_variable>())
{
}

void Solver::init(const MatX3& V, const MatX3i& F, const MatX3& V_cut, const MatX3i& F_cut, Utils::Init init, const MatX2& V_loaded, const MatX3i& F_loaded)
{
	using namespace placeholders;
	if (!full_init_done)
	{
		Utils::get_soup_and_helper_matrices(V, F, energy->separation->EVvar1, energy->separation->EVvar2, energy->separation->V2V);
		switch (init)
		{
		case Utils::Init::RANDOM:
			Utils::generate_soup_2d_random(V, F, Vs, Fs);
			break;
		case Utils::Init::ISOMETRIC:
			Utils::generate_soup_2d_iso(V, F, Vs, Fs);
			break;
		case Utils::Init::HARMONIC:
			Utils::generate_soup_2d_harmonic(V_cut, F_cut, Vs, Fs);
			break;
		case Utils::Init::LOADED:
			Vs = V_loaded;
			Fs = F_loaded;
			break;
		case Utils::Init::LOAD_FROM_FILE:
		{
			string uv_filename = igl::file_dialog_open();
			MatX3 Vst3;
			MatX2 Vst;
			MatX3i Fst;
			igl::readOBJ(uv_filename, Vst3, Fst);

			int nf = Fst.rows();
			Veci lin = Eigen::VectorXi::LinSpaced(3 * nf, 0, 3 * nf - 1);
			Fs = Eigen::Map<Mat3Xi>(lin.data(), 3, nf).transpose();

			Vst = Vst3.block(0, 0, Vst3.rows(), 2);
			Vs = MatX2(3 * Fs.rows(), 2);
			Mat32 face;
			for (int i = 0; i < Fst.rows(); ++i)
			{
				igl::slice(Vst, Fst.row(i).eval(), 1, face);
				Vs.block<3, 2>(3 * i, 0) = face;
			}
			break;
		}
		default:
			assert(false);
		}
		if (init != Utils::Init::LOADED && init != Utils::Init::LOAD_FROM_FILE)
			Utils::correct_flipped_triangles(Vs, Fs);
	}
	energy->init(F.rows(), Vs, Fs, V, F);
	m_x = Eigen::Map<Vec>(Vs.data(), Vs.rows() * Vs.cols());
	ext_x = m_x;
	eval_f = bind(&Energy::evaluate_f, energy, _1, _2);
	eval_fgh = bind(&Energy::evaluate_fgh, energy, _1, _2, _3, _4);
	internal_init();
	num_steps = 2147483647;
}

int Solver::run()
{
	is_running = true;
	halt = false;
	int steps = 0;
	do
	{
		ret = step();
		linesearch();
		update_external_mesh();
	} while ((a_parameter_was_updated || test_progress()) && !halt && ++steps < num_steps);
	is_running = false;
	cout << "solver stopped" << endl;
	return ret;
}

void Solver::stop()
{
	wait_for_param_slot();
	halt = true;
	release_param_slot();
}

void Solver::update_external_mesh()
{
	give_param_slot();
	unique_lock<shared_mutex> lock(*mesh_mutex);
	internal_update_external_mesh();
	progressed = true;
}

void Solver::get_mesh(MatX2& X)
{
	unique_lock<shared_mutex> lock(*mesh_mutex);
	Vs = Eigen::Map<MatX2>(ext_x.data(), ext_x.rows() / 2, 2);
	X = Vs;
	progressed = false;
}

void Solver::give_param_slot()
{
	a_parameter_was_updated = false;
	unique_lock<mutex> lock(*param_mutex);
	params_ready_to_update = true;
	param_cv->notify_one();
	while (wait_for_param_update)
	{
		param_cv->wait(lock);
		a_parameter_was_updated = true;
	}
	params_ready_to_update = false;
}

void Solver::wait_for_param_slot()
{
	unique_lock<mutex> lock(*param_mutex);
	wait_for_param_update = true;
	while (!params_ready_to_update && is_running)
		param_cv->wait_for(lock, chrono::milliseconds(50));
}

void Solver::release_param_slot()
{
	wait_for_param_update = false;
	param_cv->notify_one();
}