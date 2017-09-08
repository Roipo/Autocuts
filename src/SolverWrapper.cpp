#include "SolverWrapper.h"

SolverWrapper::SolverWrapper()
	:
	newton(make_shared<Newton>())
{
	solver = newton;
}

void SolverWrapper::init(const MatX3& V, const MatX3i& F, const MatX3& V_cut, const MatX3i& F_cut, Utils::Init init, const MatX2& V_loaded, const MatX3i& F_loaded)
{
	solver->init(V, F, V_cut, F_cut, init, V_loaded, F_loaded);
}

void SolverWrapper::set_bound(double new_bound)
{
	solver->wait_for_param_slot();
	solver->energy->symDirichlet->bound = new_bound;
	solver->release_param_slot();
}

void SolverWrapper::set_lambda(double new_lambda)
{
	solver->wait_for_param_slot();
	solver->energy->lambda = new_lambda;
	solver->release_param_slot();
}

void SolverWrapper::set_delta(double new_delta)
{
	solver->wait_for_param_slot();
	solver->energy->separation->delta = new_delta;
	solver->release_param_slot();
}

void SolverWrapper::set_target_triangle_position(const Mat32& pos, int hit_triangle)
{
	solver->wait_for_param_slot();
	solver->energy->position->target_pos = pos;
	solver->release_param_slot();
}

bool SolverWrapper::if_its_a_fixed_triangle_remove_and_store(int fid, Mat32& pos)
{
	solver->wait_for_param_slot();
	if (solver->energy->position->fixed_triangles.size() > 0 && solver->energy->position->fixed_triangles.find(fid) != solver->energy->position->fixed_triangles.end())
	{
		restore_constrained_face = fid;
		pos = solver->energy->position->fixed_triangles[fid];
		solver->energy->position->fixed_triangles.erase(fid);
		return true;
	}
	solver->release_param_slot();
	return false;
}

void SolverWrapper::if_moved_triangle_was_fixed_restore_constraint()
{
	if (restore_constrained_face != -1)
		add_or_remove_triangle_from_fixed_position_set(restore_constrained_face, solver->energy->position->target_pos);
	restore_constrained_face = -1;
}

void SolverWrapper::add_or_remove_triangle_from_fixed_position_set(int fid, const Mat32& pos)
{
	solver->wait_for_param_slot();
	if (solver->energy->position->fixed_triangles.size() > 0 && solver->energy->position->fixed_triangles.find(fid) != solver->energy->position->fixed_triangles.end())
		solver->energy->position->fixed_triangles.erase(fid);
	else
		solver->energy->position->fixed_triangles.insert(pair<int, Mat32>(fid, pos));
	solver->release_param_slot();
}

map<int, Mat32> SolverWrapper::add_triangles_to_position_set(const map<int, Mat32>& new_fixed_triangles, map<int, Mat32>& target_map)
{
	solver->wait_for_param_slot();
	for (auto it = new_fixed_triangles.begin(); it != new_fixed_triangles.end(); ++it)
	{
		target_map[it->first] = it->second;
	}
	map<int, Mat32> ret = target_map;
	solver->release_param_slot();
	return ret;
}

void SolverWrapper::remove_triangles_from_position_sets(const map<int, Mat32>& new_fixed_triangles)
{
	solver->wait_for_param_slot();
	for (auto it = new_fixed_triangles.begin(); it != new_fixed_triangles.end(); ++it)
	{
		solver->energy->position->fixed_triangles.erase(it->first);
		solver->energy->position->active_triangles.erase(it->first);
	}
	solver->release_param_slot();
}

void SolverWrapper::set_active_triangle(int hit_triangle)
{
	solver->wait_for_param_slot();
	solver->energy->position->active_triangle = hit_triangle;
	solver->release_param_slot();
}

void SolverWrapper::clear_fixed_triangles()
{
	solver->wait_for_param_slot();
	solver->energy->position->fixed_triangles.clear();
	solver->release_param_slot();
}

void SolverWrapper::set_position_weight(double new_pos)
{
	solver->wait_for_param_slot();
	solver->energy->pos_weight = new_pos;
	solver->release_param_slot();
}

void SolverWrapper::set_esep(const SpMat& Esept)
{
	solver->wait_for_param_slot();
	solver->energy->separation->Esept = Esept;
	solver->energy->separation->Esep = Esept.transpose();
	solver->release_param_slot();
}

void SolverWrapper::update_no_seam_constraints(int pair1, int pair2, double val)
{
	solver->wait_for_param_slot();
	solver->energy->separation->no_seam_constraints_per_pair(pair1) = val;
	solver->energy->separation->no_seam_constraints_per_pair(pair2) = val;
	solver->release_param_slot();
}

void SolverWrapper::add_indices(int idx1, int idx2, int gidx1, int gidx2, int gidx3, int gidx4)
{
	solver->wait_for_param_slot();
	solver->energy->separation->value_force_connects.push_back(idx1);
	solver->energy->separation->value_force_connects.push_back(idx2);
	solver->energy->separation->gradient_force_connects.push_back(gidx1);
	solver->energy->separation->gradient_force_connects.push_back(gidx2);
	solver->energy->separation->gradient_force_connects.push_back(gidx3);
	solver->energy->separation->gradient_force_connects.push_back(gidx4);
	solver->release_param_slot();
}

void SolverWrapper::remove_indices(int idx1, int idx2, int gidx1, int gidx2, int gidx3, int gidx4)
{
	solver->wait_for_param_slot();
	auto fit = find(solver->energy->separation->value_force_connects.begin(), solver->energy->separation->value_force_connects.end(), idx1);
	solver->energy->separation->value_force_connects.erase(fit, fit + 2); // remove idx1 and idx2
	auto git = find(solver->energy->separation->gradient_force_connects.begin(), solver->energy->separation->gradient_force_connects.end(), gidx1);
	solver->energy->separation->gradient_force_connects.erase(git, git + 4); // remove gidx1 to gidx4
	solver->release_param_slot();
}

void SolverWrapper::set_mesh_position(const MatX2& Vs_new)
{
	solver->wait_for_param_slot();
	solver->Vs = Vs_new;
	solver->m_x = Eigen::Map<const Vec>(Vs_new.data(), Vs_new.rows() * Vs_new.cols());
	solver->ext_x = Eigen::Map<const Vec>(Vs_new.data(), Vs_new.rows() * Vs_new.cols());
	solver->release_param_slot();
}

bool SolverWrapper::progressed()
{
	return solver->progressed;
}

void SolverWrapper::get_slot()
{
	solver->wait_for_param_slot();
}

void SolverWrapper::release_slot()
{
	solver->release_param_slot();
}