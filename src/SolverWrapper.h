#pragma once

#include "Solver.h"
#include "Newton.h"

class SolverWrapper
{
public:
	SolverWrapper();

	void init(const MatX3& V, const MatX3i& F, const MatX3& V_cut, const MatX3i& F_cut, Utils::Init init, const MatX2& V_loaded = MatX2(), const MatX3i& F_loaded = MatX3i());

	// Getter & setter for parameters and the mesh
	void set_bound(double new_bound);
	void set_lambda(double new_lambda);
	void set_delta(double new_delta);
	void set_target_triangle_position(const Mat32& pos, int hit_triangle);
	bool if_its_a_fixed_triangle_remove_and_store(int fid, Mat32& pos);
	void if_moved_triangle_was_fixed_restore_constraint();
	void add_or_remove_triangle_from_fixed_position_set(int fid, const Mat32& pos);
	map<int, Mat32> add_triangles_to_position_set(const map<int, Mat32>& new_fixed_triangles, map<int, Mat32>& target_map);
	void remove_triangles_from_position_sets(const map<int, Mat32>& new_fixed_triangles);
	void set_active_triangle(int hit_triangle);
	void set_position_weight(double new_pos);
	void set_esep(const SpMat& Esept);
	void update_no_seam_constraints(int pair1, int pair2, double val);

	void add_indices(int idx1, int idx2, int gidx1, int gidx2, int gidx3, int gidx4);
	void remove_indices(int idx1, int idx2, int gidx1, int gidx2, int gidx3, int gidx4);

	void set_mesh_position(const MatX2& Vs_new);
	void clear_fixed_triangles();
	bool progressed();

	void get_slot();
	void release_slot();

	// Interface to outside
	shared_ptr<Solver> solver;

private:
	// Explicit solver implementations
	shared_ptr<Newton> newton;

	// when moving a green face, restore it again to green after moving
	int restore_constrained_face = -1;
};