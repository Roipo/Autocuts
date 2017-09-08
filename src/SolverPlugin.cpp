#include "SolverPlugin.h"

#include "svg_exporter.h"

#include <queue>
#include <deque>
#include <array>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include <igl\slice.h>
#include <igl\unique.h>
#include <igl\writeOBJ.h>
#include <igl\unproject.h>
#include <igl\slice_into.h>
#include <igl\png\writePNG.h>
#include <igl\unproject_ray.h>
#include <igl\adjacency_list.h>
#include <igl\cut_mesh_simple.h>
#include <igl\file_dialog_save.h>
#include <igl\viewer\ViewerData.h>
#include <igl\ray_mesh_intersect.h>
#include <igl\all_pairs_distances.h>
#include <igl\unproject_onto_mesh.h>
#include <igl\segment_segment_intersect.h>
#include <igl\triangle_triangle_adjacency.h>
#include <igl\is_boundary_edge.h>
#include "unionFind.h"

SolverPlugin::SolverPlugin()
	: solver_wrapper(make_unique<SolverWrapper>())
{
	// get random seed for soup generation
	//srand(time(0));
	srand(0);
	init_rotation_matrices();
	// color for edge highlights
	C << 0.0, 1.0, 0.0;
	C_hover << 0.854902, 0.647059, 0.12549;
	white << 0.7, 0.7, .4;
	red << 1.0, 0.0, 0.0;
	C_merge << 51. / 255., 204. / 255., 1.0;
	zero3 << 0., 0., 0.;
	ones3 << 1., 1., 1.;
	black << 0., 0., 0.;
}

void SolverPlugin::InitPipelineStep()
{
	bar->addGroup("Mesh");
	bar->addButton("Load Scene", [&]()
	{
		string filename = igl::file_dialog_open();
		load_scene(filename);
	});
	bar->addButton("Save Scene", [&]()
	{
		string filename = igl::file_dialog_save();
		save_scene(filename);
	});
	auto uv_init_type = bar->addVariable("UV Init Type", uv_init, true);
	uv_init_type->setItems({ "RANDOM", "ISOMETRIC", "LOADED", "LOAD FROM FILE" });
	uv_init_type->setFixedWidth(140);

	bar->addGroup("Solver");

	bar->addButton("Run", [&]()
	{
		start_solver_thread();
	});
	bar->addVariable<bool>("Use num steps",
		[&](bool val)
	{
		use_num_steps = val;
		if (val)
		{
			solver_wrapper->solver->num_steps = 20;
		}
		else
		{
			solver_wrapper->solver->num_steps = INT_INF;
		}
	},
		[&]() { return use_num_steps; },
		true);
	bar->addVariable("Num steps", solver_wrapper->solver->num_steps, true)->setFixedWidth(140);
	bar->addVariable("Dist. bound", solver_wrapper->solver->energy->symDirichlet->bound)->setCallback([&](double value)
	{
		solver_wrapper->set_bound(value);
	});
	
	bar->addVariable<double>("Pos. weight",
		[&](double value) { update_energy_param(Param::POSITION_WEIGHT, value); },
		[&]() { return solver_wrapper->solver->energy->pos_weight; },
		true);

	bar->addVariable("Sep. delta", solver_wrapper->solver->energy->separation->delta)->setFixedWidth(140);
	
	add_slider();
	
	auto set_energy_type = bar->addVariable("Sep. Energy Type", solver_wrapper->solver->energy->separation->sepEType, true);
	set_energy_type->setItems({ "LOG", "QUADRATIC", "FLAT_LOG", "QUOTIENT", "QUOTIENT_NEW"});
	set_energy_type->setFixedWidth(140);
	
	bar->addGroup("Separation Measures");
	bar->addVariable("Max. Element", solver_wrapper->solver->energy->max_sep)->setFixedWidth(140);
	bar->addVariable("Gradient norm", solver_wrapper->solver->energy->grad_norm_sep)->setFixedWidth(140);

	bar->addGroup("Distortion Measures");
	bar->addVariable("Max. Element", solver_wrapper->solver->energy->max_dist)->setFixedWidth(140);
	bar->addVariable("Gradient norm", solver_wrapper->solver->energy->grad_norm_dist)->setFixedWidth(140);

	bar->addGroup("Position Measures");
	bar->addVariable("Max. Element", solver_wrapper->solver->energy->max_pos)->setFixedWidth(140);
	bar->addVariable("Gradient norm", solver_wrapper->solver->energy->grad_norm_pos)->setFixedWidth(140);

	bar->addGroup("BBox Measures");
	bar->addVariable("Max. Element", solver_wrapper->solver->energy->max_bbox)->setFixedWidth(140);
	bar->addVariable("Gradient norm", solver_wrapper->solver->energy->grad_norm_bbox)->setFixedWidth(140);

	bar->addGroup("Texture Coloring");
	bar->addVariable<nanogui::Color>("Color 1",
		[&](nanogui::Color value) { viewer->get_mesh(0).tex_col1 = value; viewer->get_mesh(0).grid_texture(); },
		[&]() { return (nanogui::Color) viewer->get_mesh(0).tex_col1; },
		true);
	bar->addVariable<nanogui::Color>("Color 2",
		[&](nanogui::Color value) { viewer->get_mesh(0).tex_col2 = value; viewer->get_mesh(0).grid_texture(); },
		[&]() { return (nanogui::Color) viewer->get_mesh(0).tex_col2; },
		true);
	bar->addVariable<bool>("Load UV from file",
		[&](bool value)
	{ 
		load_uv_from_file = value;
		if (load_uv_from_file)
		{
			uv_init = Utils::Init::LOAD_FROM_FILE;
		}
		else
		{
			uv_init = Utils::Init::ISOMETRIC;
		}
	},
		[&]() { return load_uv_from_file; },
		true)->setFixedWidth(140);
	bar_window = viewer->ngui->window();

	// Minimalistic menu (hidden in the beginning)
	mini_menu = bar->addWindow(Vec2i(10-500, 10), "Mini Menu");
	bar->addGroup("Inputs");
	bar->addVariable("Sep. delta", solver_wrapper->solver->energy->separation->delta)->setFixedWidth(140);
	add_slider();
	auto mode_items = bar->addVariable<Mode>("Mode", mode, true);
	mode_items->setItems({ "MOVE", "EDGE_CUTTING", "FACE_POSITIONING", "PAINTING", "VERTEX_CLICKING", "BBOX_DRAWING" });
	mode_items->setFixedWidth(140);
	bar->addVariable("Highlight thres", sep_thresh)->setFixedWidth(140);
	bar->addVariable("Force connect", solver_wrapper->solver->energy->separation->force_factor)->setFixedWidth(140);
	bar->addButton("Run", [&]()
	{
		start_solver_thread();
	});
	bar->addButton("Run with Result Save", [&]()
	{
		run_solver_and_store_every_iteration();
	});
	bar->addVariable("Resolution Factor", resolution_factor, true)->setFixedWidth(140);
	bar->addVariable("Max weight", max_weighting_val, true)->setFixedWidth(140);
	bar->addVariable("Step weight", weighting_step, true)->setFixedWidth(140);
	bar->addVariable("Edge weight = Average", set_edge_lenghts_to_average, true)->setFixedWidth(140);
	add_texture_slider(bar->window(), texture_size, "Texture size");
	bar->addVariable<double>("BBox weight",
		[&](double value)
		{
			solver_wrapper->get_slot();
			solver_wrapper->solver->energy->bbox_weight = value;
			solver_wrapper->release_slot();
		},
		[&]() { return solver_wrapper->solver->energy->bbox_weight; },
		true);
	bar->addVariable<bool>("Show harmonic cut",
		[&](bool value)
	{
		show_harmonic_cut = value;
		update_colors = true;
	},
		[&]() { return show_harmonic_cut; },
			true)->setFixedWidth(140);
	add_color_clamp_slider("Sep", max_sep_color_value, sep_color_clamp);
	add_color_clamp_slider("Dist", max_dist_color_value, dist_color_clamp);
	bar->addVariable("Show RGB", colorByRGB);
	bar->addVariable("Store 3d Mesh", store_3d_mesh, true)->setFixedWidth(140);
	bar->addVariable<bool>("Show separation",
		[&](bool value) { show_separation_error = value; update_colors = true; },
		[&]() { return show_separation_error; },
		true)->setFixedWidth(140);
	bar->addVariable<bool>("Show distortion",
		[&](bool value) { show_distortion_error = value; update_colors = true; },
		[&]() { return show_distortion_error; },
		true)->setFixedWidth(140);
	bar->addVariable("No-Seam power", no_seam_force, true)->setFixedWidth(140);
	bar->addVariable<bool>("Show edge highlights",
		[&](bool value) { show_fixed_3d_edge_highlights = value; update_colors = true; },
		[&]() { return show_fixed_3d_edge_highlights; },
		true)->setFixedWidth(140);
	bar->addButton("Store result", [&]() { store_result(); });
	bar->addButton("Compare to Quads", [&]() { find_new_edges_after_triangulation(); });
	bar->addVariable<bool>("Show quad mesh",
		[&](bool value)
		{
			show_quad_mesh = value;
			if (show_quad_mesh)
			{
				viewer->get_mesh(mesh_id).show_lines = false;
				viewer->core.overlay_line_width = 0.5;
			}
			else
			{
				viewer->get_mesh(mesh_id).show_lines = true;
				viewer->core.overlay_line_width = 50.;
			}
			update_colors = true;
		},
		[&]() { return show_quad_mesh; },
		true)->setFixedWidth(140);
	bar->addButton("Export UV to OBJ", [&]() { export_uv_to_obj(); });

	mini_menu_old_pos = mini_menu->position();

	viewer->screen->performLayout();

	viewer->core.is_animating = true;
	viewer->core.animation_max_fps = 30.0;
	viewer->data.object_scale = 10.0;
}

bool SolverPlugin::Load(string filename)
{
	if (solver_wrapper->solver->is_running)
		stop_solver_thread();

	bool read_obj = false;
	bool read_off = false;

	string file_format = filename.substr(filename.length() - 3, 3);
	if (file_format.compare("obj") == 0)
	{
		if (!igl::readOBJ(filename, V, F))
		{
			cerr << "Failed to load mesh: " << filename << endl;
			return false;
		}
	}
	else if (file_format.compare("off") == 0)
	{
		if (!igl::readOFF(filename, V, F))
		{
			cerr << "Failed to load mesh: " << filename << endl;
			return false;
		}
	}
	else
	{
		cerr << "Unknown file format " << filename << endl;
		return false;
	}

	init();

	mesh_loaded = true;
	
	hide_menus();
	menus_are_visible = !menus_are_visible;
	
	mesh_filename = filename;

	return true;
}

void SolverPlugin::init()
{
	if (solver_wrapper->solver->is_running)
		solver_wrapper->solver->stop();

	while (solver_wrapper->solver->is_running);

	if (V.rows() == 0 || F.rows() == 0)
		return;


	solver_wrapper->init(V, F, V_cut, F_cut, uv_init);

	if (uv_id != 0) {
		viewer->get_mesh(uv_id).clear();
		viewer->get_mesh(mesh_id).clear();
		viewer->get_mesh(vp_id).clear();
		viewer->get_mesh(xh_id).clear();
	}
	else
	{
		uv_id = viewer->add_mesh("UV");
		vp_id = viewer->add_mesh("VP");
		cs_id = viewer->add_mesh("CS");
		xh_id = viewer->add_mesh("XR");
	}

	viewer->get_mesh(uv_id).set_mesh(solver_wrapper->solver->Vs, solver_wrapper->solver->Fs);

	// setup mesh soup
	mesh_soup_F = MatX3i(F.rows(), 3);
	mesh_soup_V = MatX3(3 * F.rows(), 3);
	Mat3 face;
	RVec3i verts;
	for (int i = 0; i < F.rows(); ++i)
	{
		verts = F.row(i);
		igl::slice(V, verts, 1, face);
		mesh_soup_F.row(i) << 3 * i, 3 * i + 1, 3 * i + 2;
		mesh_soup_V.block<3, 3>(3 * i, 0) = face;
		// map from old to soup vertices
		mesh_map_to_soup[verts(0)].push_back(3 * i);
		mesh_map_to_soup[verts(1)].push_back(3 * i + 1);
		mesh_map_to_soup[verts(2)].push_back(3 * i + 2);
		// map from new vertices to "parent" vertices
		mesh_map_to_orig[3 * i].push_back(verts(0));
		mesh_map_to_orig[3 * i + 1].push_back(verts(1));
		mesh_map_to_orig[3 * i + 2].push_back(verts(2));
	}
	RGBColors = mesh_soup_V;
	RGBColors.rowwise() -= RGBColors.colwise().minCoeff();
	RGBColors *= RGBColors.colwise().maxCoeff().cwiseInverse().asDiagonal();
	//viewer->get_mesh(mesh_id).set_mesh(V, F);
	set_3d_soup_vertices_with_connected_vertices(V);

	//viewer->get_mesh(mesh_id).set_colors(MatX3::Ones(F.rows(), 3));
	viewer->get_mesh(uv_id).F;
	viewer->get_mesh(uv_id).set_colors(MatX3::Ones(solver_wrapper->solver->Fs.rows(), 3));

	// init colors of uv mesh
	uv_triangle_colors = MatX3::Ones(solver_wrapper->solver->Fs.rows(), 3);
	ones_vec = Vec::Ones(solver_wrapper->solver->Fs.rows());

	viewer->get_mesh(mesh_id).set_colors(RVec3(1., 1., 1.));

	// set adjacency
	igl::adjacency_list(F, adjacency_list);

	// set initial weights for painint
	painting_weights = MatX3::Zero(solver_wrapper->solver->Vs.rows(), 3);
	paintint_ones = MatX3::Ones(solver_wrapper->solver->Vs.rows(), 3);
	painting_colors = MatX3::Zero(solver_wrapper->solver->Vs.rows(), 3);

	Mat D;
	igl::all_pairs_distances((Mat)V, (Mat)V, true, D);
	double max_dist = D.maxCoeff();
	intersection_delta = 0.05*max_dist;

	// init edge lenghts for separation
 	find_edge_lenghts_for_separation();

	build_colored_squares();
}

void SolverPlugin::init_after_scene_load()
{
	uv_init = Utils::Init::LOADED;
	solver_wrapper->init(V, F, MatX3(), MatX3i(), uv_init, solver_wrapper->solver->Vs, solver_wrapper->solver->Fs);

	if (uv_id != 0) {
		viewer->get_mesh(uv_id).clear();
		viewer->get_mesh(mesh_id).clear();
		viewer->get_mesh(vp_id).clear();
		viewer->get_mesh(xh_id).clear();
	}
	else
	{
		uv_id = viewer->add_mesh("UV");
		vp_id = viewer->add_mesh("VP");
		cs_id = viewer->add_mesh("CS");
		xh_id = viewer->add_mesh("XR");
	}

	viewer->get_mesh(uv_id).set_mesh(solver_wrapper->solver->Vs, solver_wrapper->solver->Fs);

	// setup mesh soup
	mesh_soup_F = MatX3i(F.rows(), 3);
	mesh_soup_V = MatX3(3 * F.rows(), 3);
	Mat3 face;
	RVec3i verts;
	for (int i = 0; i < F.rows(); ++i)
	{
		verts = F.row(i);
		igl::slice(V, verts, 1, face);
		mesh_soup_F.row(i) << 3 * i, 3 * i + 1, 3 * i + 2;
		mesh_soup_V.block<3, 3>(3 * i, 0) = face;
		// map from old to soup vertices
		mesh_map_to_soup[verts(0)].push_back(3 * i);
		mesh_map_to_soup[verts(1)].push_back(3 * i + 1);
		mesh_map_to_soup[verts(2)].push_back(3 * i + 2);
		// map from new vertices to "parent" vertices
		mesh_map_to_orig[3 * i].push_back(verts(0));
		mesh_map_to_orig[3 * i + 1].push_back(verts(1));
		mesh_map_to_orig[3 * i + 2].push_back(verts(2));
	}

	//viewer->get_mesh(mesh_id).set_mesh(V, F);
	set_3d_soup_vertices_with_connected_vertices(V);

	//viewer->get_mesh(mesh_id).set_colors(MatX3::Ones(F.rows(), 3));
	viewer->get_mesh(uv_id).set_colors(MatX3::Ones(solver_wrapper->solver->Fs.rows(), 3));

	// init colors of uv mesh
	uv_triangle_colors = MatX3::Ones(solver_wrapper->solver->Fs.rows(), 3);
	ones_vec = Vec::Ones(solver_wrapper->solver->Fs.rows());

	viewer->get_mesh(mesh_id).set_colors(RVec3(1., 1., 1.));

	// set adjacency
	igl::adjacency_list(F, adjacency_list);

	// set initial weights for painint
	paintint_ones = MatX3::Ones(solver_wrapper->solver->Vs.rows(), 3);

	Mat D;
	igl::all_pairs_distances((Mat)V, (Mat)V, true, D);
	double max_dist = D.maxCoeff();
	intersection_delta = 0.05*max_dist;

	// init edge lenghts for separation
	find_edge_lenghts_for_separation();

	solver_wrapper->solver->energy->separation->update_alphas(painting_weights, max_weighting_val);

	mesh_loaded = true;

	hide_menus();
	menus_are_visible = !menus_are_visible;

	build_colored_squares();
}

void SolverPlugin::build_colored_squares()
{
	if (colored_squares_already_built)
		return;

	// place a cube at (1085, 75)
	int x = 1085;
	int y = 75;
	RVec3 cs_pos;
	viewer->core.viewport << 0, 0, 1200, 1350;
	igl::unproject(RVec3(x, viewer->screen->size()[1] - y, 0.), (viewer->core.view * viewer->get_mesh(cs_id).model).eval(), viewer->core.proj, viewer->core.viewport, cs_pos);

	// at position cs_pos place a small square
	RVec2 base_pos = cs_pos.block<1, 2>(0, 0);
	Mat82 CS_X = Mat82();
	double offset = 0.15;
	CS_X.row(0) = base_pos + RVec2(offset, offset);
	CS_X.row(1) = base_pos + RVec2(-offset, offset);
	CS_X.row(2) = base_pos + RVec2(offset, -offset);
	CS_X.row(3) = base_pos + RVec2(-offset, -offset);

	viewer->get_mesh(cs_id).add_label(base_pos - RVec2(0.65, .0), "Lambda", Vec3(0., 0., 0.));

	// place a cube at (1085, 95)
	y = 210;
	viewer->core.viewport << 0, 0, 1200, 1350;
	igl::unproject(RVec3(x, viewer->screen->size()[1] - y, 0.), (viewer->core.view * viewer->get_mesh(cs_id).model).eval(), viewer->core.proj, viewer->core.viewport, cs_pos);

	// at position cs_pos place a small square
	base_pos = cs_pos.block<1, 2>(0, 0);
	CS_X.row(4) = base_pos + RVec2(offset, offset);
	CS_X.row(5) = base_pos + RVec2(-offset, offset);
	CS_X.row(6) = base_pos + RVec2(offset, -offset);
	CS_X.row(7) = base_pos + RVec2(-offset, -offset);

	viewer->get_mesh(cs_id).add_label(base_pos - RVec2(0.65, .0), "Delta", Vec3(0., 0., 0.));

	Mat43i CS_F = Mat43i();
	CS_F << 0, 1, 2,
		1, 2, 3,
		4, 5, 6,
		5, 6, 7;

	viewer->get_mesh(cs_id).set_mesh(CS_X, CS_F);
	update_boxes_colors();
	viewer->get_mesh(cs_id).show_lines = false;

	colored_squares_already_built = true;
}

void SolverPlugin::update_boxes_colors()
{
	// first box is lambda
	double l = solver_wrapper->solver->energy->lambda;
	RVec3 lambda_col = RVec3(l, 1. - l, 0.);
	// second box is delta
	double d = pow(solver_wrapper->solver->energy->separation->delta, 0.4);
	RVec delta_col = RVec3(1. - d, d, 0.);
	// set colors
	Mat43 box_colors;
	box_colors.row(0) = lambda_col;
	box_colors.row(1) = lambda_col;
	box_colors.row(2) = delta_col;
	box_colors.row(3) = delta_col;
	// set new colors
	viewer->get_mesh(cs_id).set_colors(box_colors);
}

void SolverPlugin::find_edge_lenghts_for_separation()
{
	SpMat Esept = solver_wrapper->solver->energy->separation->Esept;
	Vec edge_lengths = Vec::Ones(Esept.cols());
	int v1, v2, num_edges, pair1, pair2;
	double edge_len;
	Vec3 v1pos, v2pos;
	double total = 0.;
	int cnt = 0;
	for (int i = 0; i < F.rows(); ++i)
	{
		Vec3i face = F.row(i);
		for (int j = 0; j < 3; ++j)
		{ // loop over all 3 triangle edges
			v1 = face(j);
			v2 = face((j + 1) % 3);
			num_edges = find_corresponding_uv_edges(v1, v2);
			if (num_edges == 2)
			{ // this is an actual separation energy pair
				v1pos = V.row(v1);
				v2pos = V.row(v2);
				edge_len = 0.5 * (v1pos - v2pos).squaredNorm();
				total += (v1pos - v2pos).norm();
				cnt++;
				pair1 = find_corresponding_pair_index(uv_edges[0].first, uv_edges[1].first);
				pair2 = find_corresponding_pair_index(uv_edges[0].second, uv_edges[1].second);
				if (pair1 == -1 || pair2 == -1)
				{
					pair1 = find_corresponding_pair_index(uv_edges[0].first, uv_edges[1].second);
					pair2 = find_corresponding_pair_index(uv_edges[0].second, uv_edges[1].first);
				}
				edge_lengths(pair1) = edge_len;
				edge_lengths(pair2) = edge_len;
			}
		}
	}
	if (set_edge_lenghts_to_average)
		edge_lengths = Vec::Constant(Esept.cols(), total / (double)cnt);
	solver_wrapper->get_slot();
	solver_wrapper->solver->energy->separation->edge_lenghts_per_pair = edge_lengths;
	solver_wrapper->release_slot();
}

void SolverPlugin::find_new_edges_after_triangulation()
{
	// 1. read quad file -> take my existing code
	string quad_name = igl::file_dialog_open();
	read_quad_mesh(quad_name);
	igl::adjacency_list(Fq, adjacency_list_quad);

	// 2. run through mesh and use adjaceny list of triangle mesh for smart lookup
	Vec3i verts;
	Mat3 face;
	// points of the edge
	RVec3 v1, v2;
	int no_seams = 0;
	for (int i = 0; i < F.rows(); ++i)
	{ // there are 4 edges to check in each quad
		verts = F.row(i);
		igl::slice(Vq, verts, 1, face);
		for (int j = 0; j < 3; ++j)
		{
			// get edge vertex indices
			ev1 = verts(j);
			ev2 = verts((j + 1) % 3);
			// check if this edge is not existing in the triangle mesh
			if (!edge_exists_in_quad_mesh(ev1, ev2) && ev1 < ev2)
			{
				// mark this edge as NO SEAM
				find_corresponding_uv_edges(ev1, ev2);
				react_to_edge_click_force_merge();
				++no_seams;
			}
		}
	}
	cout << "registered " << no_seams << " of no-seam edges" << endl;
}

void SolverPlugin::snap_soup_mat(const MatX3& Ms, const MatX3i& Fs, MatX3& Md, MatX3i& Fd)
{
	//MatX2 UV;
	//solver_wrapper->solver->get_mesh(UV);
	//MatX3i F_UV = solver_wrapper->solver->Fs;

	auto f_sep_per_pair=solver_wrapper->solver->energy->separation->f_sep_per_pair;
	auto ind2pair = solver_wrapper->solver->energy->separation->ind2pair;
	auto C2C = solver_wrapper->solver->energy->separation->C2C;
	// we have "at most" the same amount of vertices, but surely much less
	

	// a map from the old to the new (merged) vertices
	// initially every vertex is replaced by the same index
	// the upcoming loop changes this mapping, which is later
	// used as a replacement-lookup
	vector<int> fmap(Ms.rows());
	iota(fmap.begin(), fmap.end(), 0);

	int curr = 0;
	vector<bool> visited(Ms.rows(), false);
	UnionFind uf(Ms.rows());
	for (int i = 0; i < Ms.rows(); ++i)
	{
		for (SpMat::InnerIterator it(C2C, i); it; ++it) {
			int j = it.row();
			pair<int, int> p(i, j);
			if (ind2pair.count(p) == 1 && f_sep_per_pair(ind2pair.at(p)) <= 0.5)
				uf.makeSameGroup(i, j);
		}
	}
	
	std::vector< std::vector<int> > cands;
	uf.returnGroups(cands);
	Md.resize(uf.nGroups(), Eigen::NoChange);
	Md.setZero();
	for (int curr = 0; curr < cands.size(); curr++)
	{
		for (int c : cands[curr])
		{
			Md.row(curr) += Ms.row(c);
			fmap[c] = curr;
		}
		Md.row(curr) /= cands[curr].size();
	}

	// use the fmap bookmarking to update the vertex indices in F
	Fd = Fs;
	for (int i = 0; i < Fd.rows(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			Fd(i, j) = fmap[Fd(i, j)];
		}
	}
	Md.conservativeResize(uf.nGroups(), Eigen::NoChange);
}


void SolverPlugin::export_uv_to_obj()
{
	// 0.) Open file
	string filename = igl::file_dialog_save();
	ofstream out;
	out.open(filename);
	if (!out.is_open())
	{
		cout << "cant open " << filename << endl;
		return;
	}

	// 1.) snap the UV and 3d Soup
	MatX3 snap_UV, snap_3d;
	MatX3i snap_UV_F, snap_3d_F;
	// pad uv with zeros
	MatX2 UV;
	solver_wrapper->solver->get_mesh(UV);
	MatX3 padded_UV = MatX3::Zero(UV.rows(), 3);
	padded_UV.block(0, 0, UV.rows(), 2) = UV;
	snap_soup_mat(padded_UV, solver_wrapper->solver->Fs, snap_UV, snap_UV_F);
	snap_soup_mat(mesh_soup_V, mesh_soup_F, snap_3d, snap_3d_F);
	
	// 2.) write 3d coordinates v of the mesh
	for (int i = 0; i < V.rows(); ++i)
	{
		out << "v " << V.row(i) << endl;
	}

	// 3.) write the uv coordinates vt
	for (int i = 0; i < snap_UV.rows(); ++i)
	{
		out << "vt " << snap_UV.row(i) << endl;
	}

	// 4.) write the faces with f vi/vt vi/vt vi/vt,
	//     where vi and vt should match in their locations in
	//     f_3d and f_uv (index-wise (i, j))
	if (F.rows() != snap_UV_F.rows())
	{
		cout << "unclear correlation between tex and 3d vertices" << endl;
		return;
	}
	// obj is 1-indexed
// 	snap_3d_F.array() += 1;
// 	snap_UV_F.array() += 1;
	for (int i = 0; i < snap_3d_F.rows(); ++i)
	{
		RVec3i f = F.row(i).array()+1;
		RVec3i ft = snap_UV_F.row(i).array()+1;
		out << "f "
			<< f(0) << "/" << ft(0) << " "
			<< f(1) << "/" << ft(1) << " "
			<< f(2) << "/" << ft(2) << endl;
	}
	out.close();
}

void SolverPlugin::read_quad_mesh(const string& quad_name)
{
	ifstream in;
	in.open(quad_name);
	if (!in.is_open())
		return;

	int nV = 0, nF = 0;
	vector<RVec3> Vt;
	vector<RVec4i> Ft;

	string line;
	while (getline(in, line))
	{
		for (int i = 0; i < line.length(); ++i)
		{
			if (line[i] == ' ')
			{
				line[i] = ',';
			}
		}
		if (line[0] == 'v')
		{
			char* l = &line[2];
			double c1 = 0., c2 = 0., c3 = 0.;
			sscanf(l, "%lf,%lf,%lf", &c1, &c2, &c3);
			Vt.push_back(RVec3(c1, c2, c3));
		}
		else if (line[0] == 'f')
		{
			char* l = &line[2];
			int i = 2;
			string num;
			vector<int> nums;
			for (int j = 0; j < 4; ++j)
			{
				do
				{
					num.push_back(line[i]);
				} while (i < line.size() && line[++i] != ',');
				++i;
				nums.push_back(atoi(num.c_str()));
				num.clear();
			}
			Ft.push_back(RVec4i(nums[0] - 1, nums[1] - 1, nums[2] - 1, nums[3] - 1));
		}
	}

	Vq = MatX3(Vt.size(), 3);
	Fq = MatX4i(Ft.size(), 4);

	for (int i = 0; i < Vt.size(); ++i)
	{
		Vq.row(i) = Vt[i];
	}
	for (int i = 0; i < Ft.size(); ++i)
	{
		Fq.row(i) = Ft[i];
	}
}

bool SolverPlugin::edge_exists_in_quad_mesh(int p1, int p2)
{
	vector<int> neighbors = adjacency_list_quad[p1];
	for (int n : neighbors)
	{
		if (n == p2)
			return true;
	}
	return false;
}

int SolverPlugin::find_corresponding_pair_index(int i1, int i2)
{
	SpMat Esept = solver_wrapper->solver->energy->separation->Esept;
	for (int i = 0; i < Esept.outerSize(); ++i)
	{
		int idx_xi, idx_xj;
		SpMat::InnerIterator it(Esept, i);
		idx_xi = it.row();
		idx_xj = (++it).row();
		if ((idx_xi == i1 && idx_xj == i2) || (idx_xj == i1 && idx_xi == i2))
			return i;
	}
	return -1;
}

void SolverPlugin::add_texture_slider(nanogui::Window* window, double& var, const string& name)
{
	Widget *panel = new Widget(window);
	panel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 20));

	Slider* s = new Slider(panel);
	s->setValue(var);

	TextBox *textBox = new TextBox(panel);
	textBox->setValue(removeTrailingZeros(to_string(round(100 * max_texture_val * var) / 100)));
	textBox->setFixedWidth(50);
	textBox->setEditable(false);
	textBox->setFixedHeight(18);

	s->setCallback([&, textBox](float value)
	{
		textBox->setValue(removeTrailingZeros(to_string(round(100 * max_texture_val * value) / 100)));
		var = max_texture_val * value;
		update_colors = true;
	});

	var *= max_texture_val;

	viewer->ngui->addWidget(name, panel);
}

void SolverPlugin::add_slider()
{
	Widget *panel = new Widget(bar->window());
	panel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 20));

	slider = new Slider(panel);
	slider->setValue(solver_wrapper->solver->energy->lambda);
	
	TextBox *textBox = new TextBox(panel);
	textBox->setValue(removeTrailingZeros(to_string(solver_wrapper->solver->energy->lambda)));
	textBox->setFixedWidth(50);
	textBox->setEditable(true);
	textBox->setFixedHeight(18);
	
	slider->setCallback([&, textBox](float value) {
		textBox->setValue(removeTrailingZeros(to_string(round(100 * value) / 100)));
	});
	slider->setFinalCallback([&](float value) {
		update_energy_param(Param::LAMBDA, value);
	});
	
	viewer->ngui->addWidget("Lambda", panel);
}

void SolverPlugin::add_color_clamp_slider(const string& name, const shared_ptr<double>& max_value, const shared_ptr<double>& value)
{
	Widget *panel = new Widget(bar->window());
	panel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 20));

	Slider* s = new Slider(panel);
	s->setValue(*value);

	TextBox *textBox = new TextBox(panel);
	textBox->setValue(removeTrailingZeros(to_string(*value)));
	textBox->setFixedWidth(50);
	textBox->setEditable(false);
	textBox->setFixedHeight(18);

	s->setCallback([&, textBox](float new_value) {
		textBox->setValue(removeTrailingZeros(to_string(round(100 * new_value* *max_value) / 100)));
		*value = new_value * *max_value;
		update_colors = true;
	});

	viewer->ngui->addWidget(name, panel);
}

inline string SolverPlugin::removeTrailingZeros(string& s) {
	return s.erase(s.find_last_not_of('0') + 1, string::npos);
}

void SolverPlugin::start_solver_thread()
{
	cout << "start new solver" << endl;
	solver_thread = thread(&Solver::run, solver_wrapper->solver.get());
	solver_thread.detach();
}

void SolverPlugin::run_solver_and_store_every_iteration()
{
	int i = 0;
	int iters = solver_wrapper->solver->num_steps;
	solver_wrapper->solver->num_steps = 1;
	switch (uv_init)
	{
	case Utils::Init::RANDOM:
		foldername = "random";
		break;
	case Utils::Init::ISOMETRIC:
		foldername = "isometric";
		break;
	}
	string filename = "C:\\Users\\sanhuber\\Dropbox\\Paper Results\\cat_initialization\\" + foldername + "\\";
	while (i++ < iters)
	{
		solver_thread = thread(&Solver::run, solver_wrapper->solver.get());
		solver_thread.detach();
		while (!solver_wrapper->solver->progressed)
			this_thread::sleep_for(10ms);
		store_result(filename + to_string(solver_iter_counter++));
		cout << "stored result " << solver_iter_counter << endl;
		solver_wrapper->solver->progressed = false;
		update_mesh();
	}
	solver_wrapper->solver->num_steps = iters;
}

void SolverPlugin::stop_solver_thread()
{
	
}

void SolverPlugin::update_energy_param(Param p, double value, bool start_solver)
{
	switch (p)
	{
	case Param::LAMBDA:
		solver_wrapper->set_lambda(value);
		slider->setValue(value);
		slider->callback()(value);
		break;
	case Param::DELTA:
		solver_wrapper->set_delta(value);
		break;
	case Param::BOUND:
		solver_wrapper->set_bound(value);
		break;
	case Param::POSITION_WEIGHT:
		solver_wrapper->set_position_weight(value);
		break;
	default:
		assert(false && "Unknown energy parameter");
	}
	//if (!solver_wrapper->solver->is_running && start_solver)
	//	start_solver_thread();
}

void SolverPlugin::update_triangle_position(const Mat32& target_triangle_pos)
{
	solver_wrapper->set_target_triangle_position(target_triangle_pos, hit_triangle);
	if (!solver_wrapper->solver->is_running)
		start_solver_thread();
}

bool SolverPlugin::KeyDown(int key, int modifiers)
{
	switch (key)
	{
	case GLFW_KEY_LEFT:
		if (solver_wrapper->solver->energy->lambda > 0.9)
			update_energy_param(Param::LAMBDA, solver_wrapper->solver->energy->lambda - 0.01);
		else if (solver_wrapper->solver->energy->lambda >= 0.1)
			update_energy_param(Param::LAMBDA, solver_wrapper->solver->energy->lambda - 0.1);
		update_boxes_colors();
		break;
	case GLFW_KEY_RIGHT:
		if (solver_wrapper->solver->energy->lambda <= 0.98)
		{
			if (solver_wrapper->solver->energy->lambda >= 0.85)
				update_energy_param(Param::LAMBDA, solver_wrapper->solver->energy->lambda + 0.01);
			else if (solver_wrapper->solver->energy->lambda <= 0.9)
				update_energy_param(Param::LAMBDA, solver_wrapper->solver->energy->lambda + 0.1);
			update_boxes_colors();
		}
		break;
	case GLFW_KEY_DOWN:
		update_energy_param(Param::DELTA, solver_wrapper->solver->energy->separation->delta * 0.5);
		update_boxes_colors();
		break;
	case GLFW_KEY_UP:
		update_energy_param(Param::DELTA, solver_wrapper->solver->energy->separation->delta * 2.0);
		update_boxes_colors();
		break;
	case GLFW_KEY_R:
		update_energy_param(Param::BOUND, INF, false);
		update_energy_param(Param::DELTA, 1.0, false);
		update_energy_param(Param::LAMBDA, 0.0, false);
		update_energy_param(Param::POSITION_WEIGHT, 100.0, false);
		solver_wrapper->set_active_triangle(-1);
		fixed_highlighted_3d_edges.clear();
		fixed_highlighted_2d_edges.clear();
		zeroed_esep_columns.clear();
		painting_weights = MatX3::Zero(solver_wrapper->solver->Vs.rows(), 3);
		painting_colors = MatX3::Zero(solver_wrapper->solver->Vs.rows(), 3);
		solver_wrapper->clear_fixed_triangles();
		solver_iter_counter = 0;
		init();
		break;
	case GLFW_KEY_H:
		if (menus_are_visible)
			hide_menus();
		else
			show_menus();
		menus_are_visible = !menus_are_visible;
		break;
	case GLFW_KEY_J:
		if (all_menus_are_visible)
			hide_all_menus();
		else
			show_all_menus();
		all_menus_are_visible = !all_menus_are_visible;
		break;
	case GLFW_KEY_1:
	case GLFW_KEY_LEFT_ALT:
		if (mode == Mode::FACE_POSITIONING)
			reset_face_positioning();
		if (mode == Mode::EDGE_CUTTING)
			reset_edge_cutting();
		if (mode == Mode::PAINTING)
			reset_painting();
		mode = Mode::MOVE;
		break;
	case GLFW_KEY_2:
	case GLFW_KEY_TAB:
		if (mode == Mode::FACE_POSITIONING)
			reset_face_positioning();
		if (mode == Mode::PAINTING)
			reset_painting();
		mode = Mode::EDGE_CUTTING;
		break;
	case GLFW_KEY_3:
		if (mode == Mode::EDGE_CUTTING)
			reset_edge_cutting();
		if (mode == Mode::PAINTING)
			reset_painting();
		mode = Mode::FACE_POSITIONING;
		break;
	case GLFW_KEY_4:
	case GLFW_KEY_LEFT_SHIFT:
		if (mode == Mode::FACE_POSITIONING)
			reset_face_positioning();
		if (mode == Mode::EDGE_CUTTING)
			reset_edge_cutting();
		mode = Mode::PAINTING;
		find_crosshair_z_position();
		draw_crosshair();
		break;
	case GLFW_KEY_L:
		viewer->open_dialog_load_mesh();
		break;
	case GLFW_KEY_T:
		// toggle texturing
		viewer->get_mesh(mesh_id).show_texture = !viewer->get_mesh(mesh_id).show_texture;
		return true; // dont trigger fill variable
	case GLFW_KEY_I:
		// toggle intersection detection
		show_intersections = !show_intersections;
		update_colors = true;
		break;
	case GLFW_KEY_F:
		// try to fix overlapping triangles idientified in check_intersections()
		free_overlapping_triangles = !free_overlapping_triangles;
		toggle_free_overlapping_faces();
		break;
	case GLFW_KEY_B:
		if (mode == Mode::FACE_POSITIONING)
		{
			update_dot_on_mouse_pos();
			mode = Mode::BBOX_DRAWING;
		}
		else
			mode = Mode::FACE_POSITIONING;
		break;
	case GLFW_KEY_ENTER:
		solver_wrapper->solver->energy->separation->update_alphas(painting_weights, max_weighting_val);
		break;
	case GLFW_KEY_S:
		save_state_active = true;
		break;
	case GLFW_KEY_KP_1:
	case GLFW_KEY_KP_2:
	case GLFW_KEY_KP_3:
	{
		int num = key - GLFW_KEY_KP_0; // (0 = 320)
		if (save_state_active)
		{ // save the state to the mesh
			saved_states[num - 1] = V;
			state_camera_zooms[num - 1] = viewer->core.mesh_camera_zoom;
			state_normals[num - 1] = viewer->get_mesh(mesh_id).V_normals;
			save_state_active = false;
			cout << "saved camera state " << num << endl;
		}
		else
		{ // reload the mesh
			if (saved_states[num - 1].size() > 0)
			{
				V = saved_states[num - 1];
				viewer->core.mesh_camera_zoom = state_camera_zooms[num - 1];
				viewer->get_mesh(mesh_id).set_normals(state_normals[num - 1]);
				set_3d_soup_vertices_with_connected_vertices(V);
				update_colors = true;
			}
			else
			{
				cout << "no camera state " << num << endl;
			}
		}
		break;
	}
	}
	return false;
}

bool SolverPlugin::KeyUp(int key, int modifiers)
{
	switch (key)
	{
	case GLFW_KEY_LEFT_SHIFT: // painting
		reset_painting();
		mode = Mode::FACE_POSITIONING;
		break;
	case GLFW_KEY_LEFT_ALT: // moving
		mode = Mode::FACE_POSITIONING;
		break;
	case GLFW_KEY_TAB: // edge cutting
		reset_edge_cutting();
		mode = Mode::FACE_POSITIONING;
		break;
	}
	if (is_in_initial_cut_mode)
		mode = Mode::VERTEX_CLICKING;
	return false;
}

void SolverPlugin::find_crosshair_z_position()
{
	if (mouse_on_uv_side)
	{ // whole uv mesh is in z = 0 plane
		crosshair_z_pos = 1.0;
	}
	else
	{ // on 3d mesh, search for actual highest z vertex, wont change during painting
		crosshair_z_pos = viewer->get_mesh(mesh_id).V.col(2).maxCoeff() + 1.0;
	}
}

void SolverPlugin::draw_crosshair()
{
	int segments = 20;
	Vec theta = Vec::LinSpaced(segments, 0., 2. * M_PI);
	Vec3 normal(0., 0., 1);
	Mat32 nullSpace;
	nullSpace << 1, 0,
		0, 1,
		0, 0;
	MatX3 circlePoints(segments, 3);
	MatX3 circleColors(segments - 1, 3);

	RVec3 xh_pos;
	viewer->core.viewport = viewer->core.highlight_viewport;
	igl::unproject(RVec3(viewer->current_mouse_x, viewer->screen->size()[1] - viewer->current_mouse_y, 0.), (viewer->core.view * viewer->get_mesh(uv_id).model).eval(), viewer->core.proj, viewer->core.viewport, xh_pos);
	xh_center_3d << xh_pos(0), xh_pos(1), crosshair_z_pos;
	projected_xh << xh_pos(0), xh_pos(1);
	igl::repmat<RVec3, MatX3>(xh_center_3d.transpose(), segments, 1, circlePoints);

	Eigen::MatrixX3d offset(segments, 3);
	offset = Eigen::VectorXd(theta.array().cos()) * nullSpace.col(0).transpose() + Eigen::VectorXd(theta.array().sin()) * nullSpace.col(1).transpose();
	offset = offset.array().colwise() * Eigen::VectorXd::Constant(segments, xh_radius).array();
	circlePoints += offset;

	igl::repmat<Eigen::RowVector3d, Eigen::MatrixX3d>(Eigen::RowVector3d(.8, .8, .8), segments - 1, 1, circleColors);

	viewer->get_mesh(xh_id).lines.resize(0, Eigen::NoChange);
	viewer->get_mesh(xh_id).add_edges(circlePoints.topRows(segments - 1), circlePoints.bottomRows(segments - 1), circleColors);
	
	update_colors = true;
}

void SolverPlugin::reset_edge_cutting()
{
	highlighted_2d_edges.clear();
	highlighted_3d_edges.clear();
}

void SolverPlugin::reset_face_positioning()
{
	hovered_triangle = -1;
}

void SolverPlugin::reset_painting()
{
	viewer->get_mesh(xh_id).lines.resize(0, Eigen::NoChange);
}

void SolverPlugin::hide_menus()
{
	old_pos_bar = bar_window->position();
	old_viewer_bar_pos = orig_window->position();
	bar_window->setPosition(old_pos_bar - Vec2i(500, 0));
	orig_window->setPosition(old_viewer_bar_pos + Vec2i(500, 0));
	mini_menu->setPosition(mini_menu_old_pos + Vec2i(500, 0));
}

void SolverPlugin::hide_all_menus()
{
	bar_window->setPosition(bar_window->position() + Vec2i(5000, 0));
	orig_window->setPosition(orig_window->position() + Vec2i(5000, 0));
	mini_menu->setPosition(mini_menu->position() + Vec2i(5000, 0));
}

void SolverPlugin::show_menus()
{
	bar_window->setPosition(old_pos_bar);
	orig_window->setPosition(old_viewer_bar_pos);
	mini_menu->setPosition(mini_menu_old_pos);
}

void SolverPlugin::show_all_menus()
{
	bar_window->setPosition(bar_window->position() - Vec2i(5000, 0));
	orig_window->setPosition(orig_window->position() - Vec2i(5000, 0));
	mini_menu->setPosition(mini_menu->position() - Vec2i(5000, 0));
}

bool SolverPlugin::PreDraw()
{
	if (solver_wrapper->progressed() || update_colors)
		update_mesh();
	if (mouse_updated)
	{
		process_mouse_move();
		mouse_updated = false;
	}
	return false;
}

void SolverPlugin::update_mesh()
{
	// handle edge highlighting
	viewer->get_mesh(mesh_id).lines.resize(0, Eigen::NoChange);
	viewer->get_mesh(uv_id).lines.resize(0, Eigen::NoChange);
	RVec3 v1, v2;
	for (pair<int, int> vv : highlighted_3d_edges)
	{
		v1 = V.row(vv.first);
		v2 = V.row(vv.second);
		viewer->get_mesh(mesh_id).add_edges(v1, v2, C_hover);
	}
	for (pair<int, int> vv : highlighted_2d_edges)
	{
		v1 = viewer->get_mesh(uv_id).V.row(vv.first);
		v2 = viewer->get_mesh(uv_id).V.row(vv.second);
		viewer->get_mesh(uv_id).add_edges(v1, v2, C_hover);
	}
	if (show_fixed_3d_edge_highlights)
	{
		for (pair<int, int> vv : fixed_highlighted_3d_edges)
		{
			v1 = V.row(vv.first);
			v2 = V.row(vv.second);
			viewer->get_mesh(mesh_id).add_edges(v1, v2, fixed_highlight_3d_colors[vv]);
		}
		for (pair<int, int> vv : fixed_highlighted_2d_edges)
		{
			v1 = viewer->get_mesh(uv_id).V.row(vv.first);
			v2 = viewer->get_mesh(uv_id).V.row(vv.second);
			viewer->get_mesh(uv_id).add_edges(v1, v2, fixed_highlight_2d_colors[vv]);
		}
	}

	float oldwidht = viewer->core.overlay_line_width;
	viewer->core.overlay_line_width = 25.;
	if (show_quad_mesh)
	{
		// add edges as overlay for nice rendering
		int v1, v2;
		RVec4i verts;
		RVec3 c1, c2;
		for (int i = 0; i < Fq.rows(); ++i)
		{
			verts = Fq.row(i);
			for (int j = 0; j < 4; ++j)
			{
				v1 = verts(j);
				v2 = verts((j + 1) % 4);
				c1 = V.row(v1);
				c2 = V.row(v2);
				viewer->get_mesh(mesh_id).add_edges(c1, c2, black);
			}
		}
	}
	viewer->core.overlay_line_width = oldwidht;

	MatX2 newX;
	solver_wrapper->solver->get_mesh(newX);
	viewer->get_mesh(uv_id).set_mesh(newX, solver_wrapper->solver->Fs);

	// set UV of 3d mesh with newX vertices
	// prepare first for 3d mesh soup
	viewer->get_mesh(mesh_id).set_uv(texture_size * newX);

	Vec f_sep_per_edge_pair = solver_wrapper->solver->energy->separation->f_sep_per_pair;
	if (f_sep_per_edge_pair.size() == 0)
		return;

	if (colorByRGB)
		uv_triangle_colors = RGBColors;
	else
		uv_triangle_colors = MatX3::Ones(3 * F.rows(), 3);

// 	viewer->get_mesh(uv_id).set_normals(viewer->get_mesh(mesh_id).F_normals);
// 	viewer->get_mesh(uv_id).dirty |= viewer->get_mesh(uv_id).DIRTY_NORMAL;

	MatX3 uv_sep_colors = MatX3::Ones(uv_triangle_colors.rows(), 3);
	Vec vals = Vec::Zero(3 * F.rows());

	if (show_separation_error)
	{
		SpMat Esept = solver_wrapper->solver->energy->separation->Esept;
		auto pair2ind = solver_wrapper->solver->energy->separation->pair2ind;
		int idx_xi, idx_xj;
		
		for (int i = 0; i < Esept.outerSize(); ++i)
		{	
			vals(pair2ind[i].first) += f_sep_per_edge_pair(i);
			vals(pair2ind[i].second) += f_sep_per_edge_pair(i);
		}

		vals = vals / *sep_color_clamp;
		vals = 1 - vals.unaryExpr([&](double val) { return (val > 1) ? 1 : val; }).array();

		uv_sep_colors.col(1) = uv_sep_colors.col(1).cwiseProduct(vals);
		uv_sep_colors.col(2) = uv_sep_colors.col(2).cwiseProduct(vals);
	}

	MatX3 uv_dist_colors = MatX3::Ones(uv_triangle_colors.rows(), 3);
	if (show_distortion_error)
	{
		RVec dist_vals = solver_wrapper->solver->energy->symDirichlet->Efi;

		// new dist color impl
		RVec dist_err = dist_vals.transpose().array() - 4.;

		// scale to [0, dist_cutoff]
		dist_err = dist_err / *dist_color_clamp;
		// map > 1 -> 1
		dist_err= 1 - dist_err.unaryExpr([&](double val) { return (val > 1) ? 1 : val; }).array();

		// map from face to vertex coloring
		Mat3X mat_vertex_dist_vals;
		igl::repmat(dist_err, 3, 1, mat_vertex_dist_vals);
		Vec vertex_dist_vals = Eigen::Map<Vec>(mat_vertex_dist_vals.data(), 3 * mat_vertex_dist_vals.cols(), 1);
		
		uv_dist_colors.col(0) = uv_dist_colors.col(0).cwiseProduct(vertex_dist_vals);
		uv_dist_colors.col(1) = uv_dist_colors.col(1).cwiseProduct(vertex_dist_vals);
	}

#pragma omp parallel for num_threads(3)
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < uv_triangle_colors.rows(); ++j)
		{
			double val = (uv_dist_colors(j, i)+uv_sep_colors(j, i))/2;
			uv_triangle_colors(j, i) *= val;
		}
	}
	
	SpMat V2Vt = solver_wrapper->solver->energy->separation->V2Vt;
	Vec mesh_vals = Vec::Zero(V.rows());
	mesh_triangle_colors = MatX3::Ones(F.rows(), 3);
	for (int i = 0; i < V2Vt.outerSize(); ++i)
	{
		for (SpMat::InnerIterator it(V2Vt, i); it; ++it)
		{
			mesh_vals(i, 0) += vals(it.row(), 0);
		}
	}
	double max_mesh = mesh_vals.maxCoeff();
	mesh_vals /= max_mesh;

	Vec3i face;
	// draw this first, so a later activation of the hover triangle
	// also overwrites the coloring
	if (hovered_triangle != -1)
	{
		// uv
		face = solver_wrapper->solver->Fs.row(hovered_triangle);
		uv_triangle_colors.row(face(0)) << C_hover;
		uv_triangle_colors.row(face(1)) << C_hover;
		uv_triangle_colors.row(face(2)) << C_hover;
		// 3d mesh
		mesh_triangle_colors.row(hovered_triangle) << C_hover;
	}
	for (pair<int, Mat32> p : solver_wrapper->solver->energy->position->fixed_triangles)
	{
		// uv
		face = solver_wrapper->solver->Fs.row(p.first);
		uv_triangle_colors.row(face(0)) << 0, 1, 0;
		uv_triangle_colors.row(face(1)) << 0, 1, 0;
		uv_triangle_colors.row(face(2)) << 0, 1, 0;
		// 3d
		mesh_triangle_colors.row(p.first) << 0, 1, 0;
	}
	for (int p : nonOverlappingTriIDs)
	{
		// uv
		face = solver_wrapper->solver->Fs.row(p);
		uv_triangle_colors.row(face(0)) << 0, 0.91, 0.64;
		uv_triangle_colors.row(face(1)) << 0, 0.91, 0.64;
		uv_triangle_colors.row(face(2)) << 0, 0.91, 0.64;
		// 3d
		mesh_triangle_colors.row(p) << 0, 0.91, 0.64;
	}
	if (hit_triangle != -1)
	{
		// uv
		face = solver_wrapper->solver->Fs.row(hit_triangle);
		uv_triangle_colors.row(face(0)) << 0, 0, 1;
		uv_triangle_colors.row(face(1)) << 0, 0, 1;
		uv_triangle_colors.row(face(2)) << 0, 0, 1;
		// 3d
		mesh_triangle_colors.row(hit_triangle) << 0, 0, 1;
	}

	viewer->get_mesh(mesh_id).points.resize(0, Eigen::NoChange);
	
	if (show_harmonic_cut)
	{
		if (initial_cut.first != -1)
		{
			RVec3 pos = viewer->get_mesh(mesh_id).V.row(initial_cut.first);
			viewer->get_mesh(mesh_id).add_points((Mat)pos, (Mat)red);
		}
		if (initial_cut.second != -1)
		{
			RVec3 pos = viewer->get_mesh(mesh_id).V.row(initial_cut.second);
			viewer->get_mesh(mesh_id).add_points((Mat)pos, (Mat)red);
		}
		if (hovered_vertex != -1)
		{
			RVec3 pos = viewer->get_mesh(mesh_id).V.row(hovered_vertex);
			viewer->get_mesh(mesh_id).add_points((Mat)pos, (Mat)red);
		}
		if (path.size() != 0)
		{
			for (int i = 1; i < path.size(); ++i)
			{
				v1 = viewer->get_mesh(mesh_id).V.row(mesh_map_to_soup[path[i - 1]][0]);
				v2 = viewer->get_mesh(mesh_id).V.row(mesh_map_to_soup[path[i]][0]);
				viewer->get_mesh(mesh_id).add_edges(v1, v2, red);
			}
		}
	}

	if (show_intersections)
	{
		// TODO: Above calculations are not needed when this is enabled
		check_intersections();
		uv_triangle_colors = MatX3::Zero(3 * F.rows(), 3);
		// dont use a loop here!
		Vec3 o3;
		o3 << 1., 1., 1.;
		for (int i = 0; i < overlapping_triangles.rows(); ++i)
		{
			if (overlapping_triangles(i))
			{
				uv_triangle_colors.block<3, 1>(3 * i, 0) = o3;
			}
		}
	}

	viewer->get_mesh(uv_id).points.resize(0, Eigen::NoChange);
	RVec3 bb_color = mode == Mode::BBOX_DRAWING ? red : C_hover;
	if (mode == Mode::BBOX_DRAWING && mouse_on_uv_side && !bbox_exists)
	{
		viewer->get_mesh(uv_id).add_points((Mat)dot_on_mouse_pos, (Mat)bb_color);
	}

	for (int i = 0; i < bbox_corners.size(); ++i)
	{
		viewer->get_mesh(uv_id).add_points((Mat)bbox_corners[i].first, (Mat)bb_color);

		if (bbox_corners[i].second != RVec2(-1., -1.))
		{ // both corners found, draw full rectangle
			viewer->get_mesh(uv_id).add_points((Mat)bbox_corners[i].second, (Mat)bb_color);

			viewer->get_mesh(uv_id).add_edges(top_left[i], bottom_left[i], bb_color);
			viewer->get_mesh(uv_id).add_edges(bottom_left[i], bottom_right[i], bb_color);
			viewer->get_mesh(uv_id).add_edges(bottom_right[i], top_right[i], bb_color);
			viewer->get_mesh(uv_id).add_edges(top_right[i], top_left[i], bb_color);
		}
	}

	if (mode == Mode::PAINTING)
	{
		MatX3 colors = painting_colors.array() / max_weighting_val;
		viewer->get_mesh(uv_id).set_colors(paintint_ones - colors);
		viewer->get_mesh(mesh_id).set_colors(paintint_ones - colors);		
	}
	else
	{
		viewer->get_mesh(uv_id).set_colors(uv_triangle_colors);
		viewer->get_mesh(mesh_id).set_colors(uv_triangle_colors);
	}

	viewer->get_mesh(uv_id).dirty |= viewer->get_mesh(uv_id).DIRTY_AMBIENT;
	viewer->get_mesh(mesh_id).dirty |= viewer->get_mesh(mesh_id).DIRTY_AMBIENT;

	update_colors = false;
}

void SolverPlugin::map_to_soup_matrix(const MatX2 in, const MatX3i& F, MatX2& out)
{
	Mat32 face;
	RVec3i verts;
	out = MatX2(3 * F.rows(), 2);
	for (int i = 0; i < F.rows(); ++i)
	{
		verts = F.row(i);
		igl::slice(in, verts, 1, face);
		out.block<3, 2>(3 * i, 0) = face;
	}
}

void SolverPlugin::check_intersections()
{
	MatX3 Vs = viewer->get_mesh(uv_id).V;
	MatX3i Fs = viewer->get_mesh(uv_id).F;
	Mati intersections = Mati::Zero(Fs.rows(), Fs.rows());

	for (int i = 0; i < Fs.rows(); ++i)
	{
		Edge e1_1, e1_2, e1_3;
		find_edges_of_triangle(i, Vs, Fs, e1_1, e1_2, e1_3);

		for (int j = 0; j < Fs.rows(); ++j)
		{
			Edge e2_1, e2_2, e2_3;
			find_edges_of_triangle(j, Vs, Fs, e2_1, e2_2, e2_3);

			bool b1 = test_intersection(e1_1, e2_1) || test_intersection(e1_1, e2_2) || test_intersection(e1_1, e2_3);
			bool b2 = test_intersection(e1_2, e2_1) || test_intersection(e1_2, e2_2) || test_intersection(e1_2, e2_3);
			bool b3 = test_intersection(e1_3, e2_1) || test_intersection(e1_3, e2_2) || test_intersection(e1_3, e2_3);

			intersections(i, j) = b1 || b2 || b3;
		}
	}
	Mati without_self_intersection = intersections - Mati::Identity(Fs.rows(), Fs.rows());
	Veci without_self_ints_sum = without_self_intersection.rowwise().sum();
	overlapping_triangles = without_self_ints_sum.array() > 0;
}

bool SolverPlugin::test_intersection(const Edge& e1, const Edge& e2)
{
	double t, u, eps = 1e-3;
	bool intersect = igl::segments_intersect(e1.first, e1.second, e2.first, e2.second, t, u, eps);
	double diff = 0.05;
	double lo = diff, hi = 1. - diff; // dont count intersections when fitting (almost) perfectly
	return intersect && t >= lo && t <= hi && u >= lo && u <= hi;
}

void SolverPlugin::toggle_free_overlapping_faces()
{
	if (overlapping_triangles.sum() == 0)
		return;

	// run through all faces, if one is overlapping another, free it
	if (free_overlapping_triangles)
	{
		MatX3i Fs = viewer->get_mesh(uv_id).F;
		for (int i = 0; i < overlapping_triangles.rows(); ++i)
		{
			if (overlapping_triangles(i))
			{
				toggle_free_face(Fs.row(i).eval());
			}
		}
	}
	else
	{ // recollect freed triangles
		while (fixed_highlighted_3d_edges.size() > 0)
		{
			// choose to delete this edge
			pair<int, int> head = fixed_highlighted_3d_edges.front();
			ev1 = head.first;
			ev2 = head.second;

			num_uv_edges_found = find_corresponding_uv_edges(ev1, ev2);

			react_to_edge_click();
		}
	}
}

void SolverPlugin::toggle_free_face(const RVec3i& face)
{
	// free face by freeing all its 3 edges
	int uv_ev1, uv_ev2;
	for (int i = 0; i < 3; ++i)
	{
		uv_ev1 = face(i);
		uv_ev2 = face((i + 1) % 3);

		find_corresponding_mesh_edge(uv_ev1, uv_ev2, ev1, ev2);
		num_uv_edges_found = find_corresponding_uv_edges(ev1, ev2);

		react_to_edge_click();
	}
}

void SolverPlugin::find_edges_of_triangle(int face, const MatX3& Vs, const MatX3i& Fs, Edge& e1, Edge& e2, Edge& e3)
{
	// get triangle t
	Mat3 t;
	igl::slice(Vs, Fs.row(face).eval(), 1, t);

	// find the corners
	Vec3 v1 = t.row(0);
	Vec3 v2 = t.row(1);
	Vec3 v3 = t.row(2);

	// find the edges
	e1 = Edge(v1, v2 - v1);
	e2 = Edge(v2, v3 - v2);
	e3 = Edge(v3, v1 - v3);
}

bool SolverPlugin::line_intersects_triangle(const Edge& e, const Mat3& T)
{
	// e.first  = origin
	// e.second = direction

	// Define corners
	Vec3 v1 = T.row(0), v2 = T.row(1), v3 = T.row(2);
	// find three edges of triangle
	Edge t1 = Edge(v1, v2 - v1), t2 = Edge(v2, v3 - v2), t3 = Edge(v3, v1 - v3);

	double tt1, tt2, tt3, uu1, uu2, uu3, eps = 1e-6;
	bool b1 = igl::segments_intersect(e.first, e.second, t1.first, t1.second, tt1, uu1, eps);
	bool b2 = igl::segments_intersect(e.first, e.second, t2.first, t3.second, tt2, uu2, eps);
	bool b3 = igl::segments_intersect(e.first, e.second, t2.first, t3.second, tt3, uu3, eps);
	bool f1 = tt1 <= 1. && uu1 <= 1.;
	bool f2 = tt2 <= 1. && uu2 <= 1.;
	bool f3 = tt3 <= 1. && uu3 <= 1.;
	return (b1 && f1) || (b2 && f2) || (b3 && f3);
}

bool SolverPlugin::MouseDown(int button, int modifier)
{
	update_highlights = false;
	viewer->core.viewport << 0, 0, 2400, 1350;
	igl::unproject(RVec3(viewer->down_mouse_x, viewer->screen->size()[1] - viewer->down_mouse_y, 0.), (viewer->core.view * viewer->get_mesh(uv_id).model).eval(), viewer->core.proj, viewer->core.viewport, projected_mouse_down);
	mesh_pos_down = V;
	igl::viewer::ViewerData dt = viewer->get_mesh(uv_id);
	uv_mesh_pos_down = viewer->get_mesh(uv_id).V;
	mesh_3d_normals_down = viewer->get_mesh(0).V_normals;

	bool LEFT_ALT = modifier == 4;
	bool LEFT_CTRL = modifier == 2;
	MOUSE_LEFT = button == 0;
	MOUSE_MID = button == 1;
	MOUSE_RIGHT = button == 2;

	switch (mode)
	{
	case Mode::MOVE:
		if (mouse_on_uv_side)
		{
			if (MOUSE_MID)
			{
				// if the solver is running, stop it during translation, and release it when releaseing the mouse
				if (solver_wrapper->solver->is_running)
				{ // stop and remember to release when done
					solver_wrapper->get_slot();
					release_when_translation_done = true;
				}
				uv_translation_enabled = true;
			}
		}
		else
		{ // mouse on 3d mesh side
			if (MOUSE_LEFT)
			{
				rotation_enabled = true;
			}
			if (MOUSE_MID)
			{
				translation_enabled = true;
			}
		}
		break;
	case Mode::EDGE_CUTTING:
		if (draws_highlighted_edge)
		{ // if an edge is highlihgted, make it clickable
			if (MOUSE_LEFT)
			{
				react_to_edge_click();
			}
			if (MOUSE_RIGHT)
			{
				react_to_edge_click_force_merge();
			}
		}
		break;
	case Mode::FACE_POSITIONING:
		if (mouse_on_uv_side)
		{
			if (mouse_over_2d_mesh)
			{
				if (MOUSE_LEFT)
				{
					if (LEFT_CTRL)
					{
						add_triangle_to_fixed_position_set();
					}	
					else
					{
						set_hit_triangle_as_active();
						move_triangle_enabled = true;
					}
				}
				else if (MOUSE_RIGHT)
				{
					largest_overlapping_part_enabled = true;
					find_largest_non_overlapping_set();
					freePatchBoundaryEdges();
				}
			}
		}
		break;
	case Mode::PAINTING:
		if (mouse_on_uv_side)
		{
			if (MOUSE_LEFT)
			{
				new_face_hits_fixed.clear();
				new_face_hits_remove.clear();
				new_face_hits_active.clear();
			}
		}
		applying_weight_enabled = true;
		break;
	case Mode::BBOX_DRAWING:
	{
		int ibox, icorner;
		bool intersects = check_if_intersects_with_a_box_corner(ibox, icorner);
		if (intersects)
		{ // we did hit a box dot, so we want to move that
			translate_box_nr = ibox;
			if (icorner == 1)
				translate_bbox_corner_1 = true;
			else
				translate_bbox_corner_2 = true;
		}
		else
		{ // there was no intersection with a dot, so we add a new box to the set
			// this is done by adding two corners consecutively
			if (is_adding_first_corner)
			{ // we add first corner
				bbox_corners.push_back(pair<RVec2, RVec2>(RVec2(-1., -1.), RVec2(-1., -1.)));
				bbox_corners.back().first << dot_on_mouse_pos(0), dot_on_mouse_pos(1);
				is_adding_first_corner = false;
			}
			else
			{ // we add the second corner
				bbox_corners.back().second << dot_on_mouse_pos(0), dot_on_mouse_pos(1);
				down_corners.push_back(bbox_corners.back());
				add_box();
				is_adding_first_corner = true;
			}
		}
		update_colors = true;
		break;
	}
	default:
		assert(false && "Unknown mode");
	}
	return true;
}

bool SolverPlugin::MouseUp(int button, int modifier)
{
	if (rotation_enabled || translation_enabled || move_triangle_enabled || hit_triangle != -1)
		update_colors = true;
	rotation_enabled = false;
	translation_enabled = false;
	if (move_triangle_enabled)
	{
		solver_wrapper->if_moved_triangle_was_fixed_restore_constraint();
	}
	move_triangle_enabled = false;
	applying_weight_enabled = false;
	update_highlights = true;
	if (uv_translation_enabled)
	{
		MatX2 Vs_new = viewer->get_mesh(uv_id).V.block(0, 0, solver_wrapper->solver->Vs.rows(), 2);
		solver_wrapper->solver->Vs = Vs_new;
		solver_wrapper->solver->m_x = Eigen::Map<const Vec>(Vs_new.data(), Vs_new.rows() * Vs_new.cols());
		solver_wrapper->solver->ext_x = Eigen::Map<const Vec>(Vs_new.data(), Vs_new.rows() * Vs_new.cols());
		if (release_when_translation_done)
			solver_wrapper->release_slot();
	}
	if (translate_bbox_corner_1 || translate_bbox_corner_2)
	{
		down_corners[translate_box_nr] = bbox_corners[translate_box_nr];
		translate_box_nr = -1;
		translate_bbox_corner_1 = false;
		translate_bbox_corner_2 = false;
	}
	uv_translation_enabled = false;
	if (button == 2 && largest_overlapping_part_enabled)
	{
		nonOverlappingTriIDs.clear();
		freePatchBoundaryEdges();
		largest_overlapping_part_enabled = false;
	}

	hit_triangle = -1;
	solver_wrapper->set_active_triangle(hit_triangle);
	return true;
}

bool SolverPlugin::MouseScroll(float delta_y)
{
	if (mode == Mode::PAINTING)
	{
		xh_radius *= delta_y > 0 ? 0.9 : 1.1;
	}
	else
	{
		if (mouse_on_uv_side)
		{ // on uv side
			viewer->core.uv_camera_zoom += delta_y > 0 ? .1 : -.1;
		}
		else
		{ // on mesh side
			viewer->core.mesh_camera_zoom += delta_y > 0 ? .1 : -.1;
		}
	}
	return true;
}

bool SolverPlugin::MouseMove(int mouse_x, int mouse_y)
{
	curr_mouse_x = mouse_x;
	curr_mouse_y = mouse_y;
	mouse_updated = true;
	return true;
}

void SolverPlugin::process_mouse_move()
{
	int mouse_x = curr_mouse_x;
	int mouse_y = curr_mouse_y;
	add_viewport_highlights(mouse_x);
	RVec3 curr_mouse_pos, curr_screen_space_mouse_pos;
	curr_screen_space_mouse_pos << mouse_x, viewer->screen->size()[1] - mouse_y, 0.;
	viewer->core.viewport << 0, 0, 2400, 1350;
	igl::unproject(curr_screen_space_mouse_pos, (viewer->core.view * viewer->get_mesh(uv_id).model).eval(), viewer->core.proj, viewer->core.viewport, curr_mouse_pos);
	double x_diff = 2 * (curr_mouse_pos(0) - projected_mouse_down(0));
	double y_diff = curr_mouse_pos(1) - projected_mouse_down(1);
	if (rotation_enabled)
	{
		rotate(-y_diff, -x_diff);
		return;
	}
	if (translation_enabled)
	{
		translate(x_diff, y_diff);
		return;
	}
	if (move_triangle_enabled)
	{
		translate_triangle(x_diff, y_diff);
		update_triangle_position(moving_triangle_pos);
		return;
	}
	if (uv_translation_enabled)
	{
		translate_uv_mesh(x_diff, y_diff);
		return;
	}
	if (uv_translation_enabled)
	{
		translate_uv_mesh(x_diff, y_diff);
		return;
	}
	// test if we are close to an edge of the 3d mesh
	if (mode == Mode::EDGE_CUTTING)
	{
		if (close_to_3d_mesh_edge() || close_to_2d_mesh_edge())
			return;
	}
	if (mode == Mode::FACE_POSITIONING)
	{
		if (mouse_on_uv_side && close_to_2d_face())
		{
			return;
		}
		else if (!mouse_on_uv_side && close_to_3d_face())
		{
			return;
		}
	}
	if (mode == Mode::PAINTING)
	{
		draw_crosshair();
		if (applying_weight_enabled)
		{
			apply_weight_onto_hit_corners();
		}
		return;
	}
	if (mode == Mode::VERTEX_CLICKING)
	{
		if (!mouse_on_uv_side && close_to_3d_vertex())
		{
			return;
		}
	}
	if (mode == Mode::BBOX_DRAWING)
	{
		if (mouse_on_uv_side)
		{
			if (translate_box_nr == -1)
			{ // no box corner is selected, just draw and be ready for a mouse down
				draw_dot_on_mouse_location();
				return;
			}
			else
			{ // we did hit a box corner and want to move it
				if (translate_bbox_corner_1)
				{
					bbox_corners[translate_box_nr].first = down_corners[translate_box_nr].first + RVec2(x_diff, y_diff);
				}
				else if (translate_bbox_corner_2)
				{
					bbox_corners[translate_box_nr].second = down_corners[translate_box_nr].second + RVec2(x_diff, y_diff);
				}
				if (translate_bbox_corner_1 || translate_bbox_corner_2)
				{
					update_box(translate_box_nr);
					update_colors = true;
					return;
				}
			}
		}
	}
	return;
}

void SolverPlugin::find_boundary_on_3d_mesh()
{
	// find parents of the viewer soup mesh in the original mesh
	pair<int, int> orig_verts = pair<int, int>(mesh_map_to_orig[initial_cut.first][0], mesh_map_to_orig[initial_cut.second][0]);
	int src = orig_verts.first;
	int dst = orig_verts.second;
	
	// breadth-first-search
	int n = V.rows();
	vector<bool> visited(n, false);
	visited[src] = true;
	vector<int> pi(n, -1);
	vector<int> d(n, 0);
	d[src] = 0;

	deque<int> Q;
	Q.push_back(src);

	while (Q.size() != 0)
	{
		int u = Q.front();
		Q.pop_front();

		for (int v : adjacency_list[u])
		{ // all neighbors
			if (!visited[v])
			{
				visited[v] = true;
				d[v] = d[u] + 1;
				pi[v] = u;
				Q.push_back(v);
			}
		}
		visited[u] = true;
	}

	// go back from dst to src
	int curr = dst;
	while (curr != src)
	{
		path.push_back(curr);
		curr = pi[curr];
	}
	path.push_back(src);
	int i = path.size() - 2;
	vector<int> back_path;
	while (i > 0)
	{
		back_path.push_back(path[i--]);
	}
	path.insert(path.end(), back_path.begin(), back_path.end());

	V_cut = V;
	F_cut = F;
	igl::cut_mesh(V_cut, F_cut, path);
}

void SolverPlugin::apply_weight_onto_hit_corners()
{
	// find all uv vertices, that lie in the xrosshair circle
	MatX3 UV;
	if (mouse_on_uv_side)
		UV = viewer->get_mesh(uv_id).V;
	else
		UV = mesh_soup_V;

	Vec3 origin = Vec3(0., 0., 20.);
	int threads = omp_get_max_threads();
	vector<list<pair<int, double>>> hits = vector<list<pair<int, double>>>(threads, list<pair<int, double>>());
	Vec2 c = projected_xh;
#pragma omp parallel for num_threads(threads)
	for (int i = 0; i < UV.rows(); ++i)
	{
		int tid = omp_get_thread_num();
		Vec3 q = UV.row(i);
		if ((q(0) - c(0))*(q(0) - c(0)) + (q(1) - c(1))*(q(1) - c(1)) - xh_radius*xh_radius < 0)
		{
			hits[tid].push_back(pair<int, double>(i, (q - origin).norm()));
		}
	}
	list<pair<int, double>> allhits;
	for (int i = 0; i < threads; ++i)
		allhits.splice(allhits.end(), hits[i]);

	if (allhits.size() == 0)
		return;

	list<int> filtered_hits;
	vector<int> vertices;
	// filter out first hits, if were on 3d
	if (!mouse_on_uv_side)
	{		
		MatX3 normals = viewer->get_mesh(mesh_id).V_normals;
		RVec3 xhc3dt = xh_center_3d.transpose();
		for (auto it = allhits.begin(); it != allhits.end(); ++it)
		{
			RVec3 n = normals.row(it->first);
			RVec uv = (xhc3dt - UV.row(it->first)).normalized();
			if (n.dot(uv) > 0.5)
				vertices.push_back(it->first);
		}
	}
	else
	{
		for (auto it = allhits.begin(); it != allhits.end(); ++it)
			vertices.push_back(it->first);
	}

	if (MOUSE_MID)
	{ // erase
		for (int i = 0; i < vertices.size(); ++i)
		{
			painting_weights.row(vertices[i]) = zero3;
			painting_colors.row(vertices[i]) = zero3;
			painting_colors.row(vertices[i]) = zero3;
		}
	}
	else
	{
		int col = MOUSE_LEFT ? 0 : 1;
		Vec2i color_cols;
		if (col == 0)
			color_cols << 1, 2;
		else if (col == 1)
			color_cols << 0, 1;
		int othercol = col == 0 ? 1 : 0;
		for (int i = 0; i < vertices.size(); ++i)
		{
			if (painting_weights(vertices[i], col) < max_weighting_val)
			{
				painting_weights(vertices[i], col) += weighting_step;
				painting_colors(vertices[i], col) += weighting_step;
				// no color/weight mixing allowed
				painting_weights(vertices[i], othercol) = 0.;
				painting_colors(vertices[i], othercol) = 0.;
			}
		}
	}

	update_colors = true;
}

int SolverPlugin::find_corresponding_uv_edges(int v1, int v2)
{
	SpMat V2Vt = solver_wrapper->solver->energy->separation->V2Vt;
	vector<int> v1_cand, v2_cand;
	// select candidates
	for (int i = 0; i < V2Vt.outerSize(); ++i)
	{
		SpMat::InnerIterator it(V2Vt, i);
		if (it.col() == v1)
		{
			for (; it; ++it)
				v1_cand.push_back(it.row());
			if (v2_cand.size() > 0)
				break; // we found both
		}
		else if (it.col() == v2)
		{
			for (; it; ++it)
				v2_cand.push_back(it.row());
			if (v1_cand.size() > 0)
				break; // we found both
		}
		continue;
	}
	// find <= 2 existing edges
	uv_edges.clear();
	for (int v1c : v1_cand)
	{
		for (int v2c : v2_cand)
		{
			if (floor(v1c / 3.0) == floor(v2c / 3.0))
				v1c < v2c ? uv_edges.push_back(pair<int, int>(v1c, v2c)) : uv_edges.push_back(pair<int, int>(v2c, v1c));
		}
	}
	return uv_edges.size();
}

void SolverPlugin::find_corresponding_mesh_edge(int uv_ev1, int uv_ev2, int& ev1, int& ev2)
{
	SpMat V2V = solver_wrapper->solver->energy->separation->V2V;
	SpMat::InnerIterator it(V2V, uv_ev1);
	ev1 = it.row();
	SpMat::InnerIterator it2(V2V, uv_ev2);
	ev2 = it2.row();
}

void SolverPlugin::react_to_edge_click()
{
	// be sure to find the first hit asap
	bool check_both = uv_edges.size() == 2;
	if (check_both && uv_edges[0].first > uv_edges[1].first)
		swap(uv_edges[0], uv_edges[1]);

	if (remove(fixed_highlighted_3d_edges, pair<int, int>(ev1, ev2)))
	{
		fixed_highlighted_3d_edges.remove(pair<int, int>(ev1, ev2));

		remove(fixed_highlighted_2d_edges, pair<int, int>(uv_edges[0].first, uv_edges[0].second));
		fixed_highlighted_2d_edges.remove(pair<int, int>(uv_edges[0].first, uv_edges[0].second));
		if (uv_edges.size() == 2)
		{
			remove(fixed_highlighted_2d_edges, pair<int, int>(uv_edges[1].first, uv_edges[1].second));
			fixed_highlighted_2d_edges.remove(pair<int, int>(uv_edges[1].first, uv_edges[1].second));
		}
	}
	else
	{
		fixed_highlighted_3d_edges.push_back(pair<int, int>(ev1, ev2));
		fixed_highlight_3d_colors[pair<int, int>(ev1, ev2)] = C;

		fixed_highlighted_2d_edges.push_back(pair<int, int>(uv_edges[0].first, uv_edges[0].second));
		fixed_highlight_2d_colors[pair<int, int>(uv_edges[0].first, uv_edges[0].second)] = C;
		if (uv_edges.size() == 2)
		{
			fixed_highlighted_2d_edges.push_back(pair<int, int>(uv_edges[1].first, uv_edges[1].second));
			fixed_highlight_2d_colors[pair<int, int>(uv_edges[1].first, uv_edges[1].second)] = C;
		}
	}
	update_colors = true;

	if (uv_edges.size() == 1)
		return;

	Vec & disconnect_alphas = solver_wrapper->solver->energy->separation->disconnect_alphas;
	auto & const ind2pair = solver_wrapper->solver->energy->separation->ind2pair;
	// each 3d edge is split into 2 uv edges, so finding one is enough
	if ((zeroed_map_it = zeroed_esep_columns.find(pair<int, int>(ev1, ev2))) != zeroed_esep_columns.end())
	{ // the edge was already clicked -> deactivate (dont ignore separation anymore, reinsert values 1 & -1
		// fill first one always
		pair<int, int> cols = zeroed_map_it->second;
// 		SpMat::InnerIterator it1(Esept, cols.first);
// 		it1.valueRef() = 1;
//		(++it1).valueRef() = -1;
		disconnect_alphas(cols.first) = 1;
		disconnect_alphas(cols.second) = 1;
		zeroed_esep_columns.erase(pair<int, int>(ev1, ev2));
	}
	else
	{ // the edge was not active -> add it to the zeroed map
		pair<int, int> p1(uv_edges[0].first, uv_edges[1].first);
		pair<int, int> p2(uv_edges[0].second, uv_edges[1].second);
		if (ind2pair.count(p1) == 0)
		{
			p1 = make_pair(uv_edges[0].first, uv_edges[0].second);
			p2 = make_pair(uv_edges[1].first, uv_edges[1].second);
			if (ind2pair.count(p1) == 0)
			{
				p1 = make_pair(uv_edges[0].first, uv_edges[1].second);
				p2 = make_pair(uv_edges[1].first, uv_edges[0].second);
			}
		}

		disconnect_alphas(ind2pair[p1]) = 0;
		disconnect_alphas(ind2pair[p2]) = 0;
		zeroed_esep_columns[pair<int, int>(ev1, ev2)] = make_pair(ind2pair[p1], ind2pair[p2]);
	}
}

void SolverPlugin::react_to_edge_click_force_merge()
{
	if (uv_edges.size() == 1)
		return;

	// be sure to find the first hit asap
	bool check_both = uv_edges.size() == 2;
	if (check_both && uv_edges[0].first > uv_edges[1].first)
		swap(uv_edges[0], uv_edges[1]);

	bool remove_edges = false;
	if (remove(fixed_highlighted_3d_edges, pair<int, int>(ev1, ev2)))
	{
		remove_edges = true;
		fixed_highlighted_3d_edges.remove(pair<int, int>(ev1, ev2));

		remove(fixed_highlighted_2d_edges, pair<int, int>(uv_edges[0].first, uv_edges[0].second));
		fixed_highlighted_2d_edges.remove(pair<int, int>(uv_edges[0].first, uv_edges[0].second));
		if (uv_edges.size() == 2)
		{
			remove(fixed_highlighted_2d_edges, pair<int, int>(uv_edges[1].first, uv_edges[1].second));
			fixed_highlighted_2d_edges.remove(pair<int, int>(uv_edges[1].first, uv_edges[1].second));
		}
	}
	else
	{
		fixed_highlighted_3d_edges.push_back(pair<int, int>(ev1, ev2));
		fixed_highlight_3d_colors[pair<int, int>(ev1, ev2)] = C_merge;

		fixed_highlighted_2d_edges.push_back(pair<int, int>(uv_edges[0].first, uv_edges[0].second));
		fixed_highlight_2d_colors[pair<int, int>(uv_edges[0].first, uv_edges[0].second)] = C_merge;
		if (uv_edges.size() == 2)
		{
			fixed_highlighted_2d_edges.push_back(pair<int, int>(uv_edges[1].first, uv_edges[1].second));
			fixed_highlight_2d_colors[pair<int, int>(uv_edges[1].first, uv_edges[1].second)] = C_merge;
		}
	}
	update_colors = true;

	// go through Esept and find the two pair indices of the clicked edges to give these pairs a higher weight
	SpMat Esept = solver_wrapper->solver->energy->separation->Esept;
	int idx_xi, idx_xj, pair1 = -1, pair2 = -1;
	for (int i = 0; i < Esept.outerSize(); ++i)
	{
		SpMat::InnerIterator it(Esept, i);
		idx_xi = it.row();
		idx_xj = (++it).row();
		if (idx_xi == uv_edges[0].first || idx_xi == uv_edges[0].second &&
			idx_xj == uv_edges[1].first || idx_xj == uv_edges[1].second)
		{
			if (pair1 == -1)
			{
				pair1 = i;
				if (pair2 != -1)
					break;
			}
			else if (pair2 == -1)
			{
				pair2 = i;
				if (pair1 != -1)
					break;
			}
		}
	}
	solver_wrapper->update_no_seam_constraints(pair1, pair2, remove_edges ? 0. : no_seam_force);
}

bool SolverPlugin::close_to_3d_mesh_edge()
{
	int fid;
	Eigen::Vector3f bc;
	viewer->core.camera_zoom = 1.0;
	viewer->core.viewport << 1200, 0, 1200, 1350;
	double x = viewer->current_mouse_x;
	double y = viewer->core.viewport(3) - viewer->current_mouse_y;
	double thresh = 0.25;
	draws_highlighted_edge = false;
	viewer->get_mesh(mesh_id).dirty |= viewer->get_mesh(mesh_id).DIRTY_OVERLAY_LINES;
	if (igl::unproject_onto_mesh(Vec2f(x, y), viewer->core.view * viewer->get_mesh(mesh_id).model,
		viewer->core.proj, viewer->core.viewport, V, F, fid, bc))
	{
		mouse_over_3d_mesh = true;
		if (bc(0) < thresh)
		{
			ev1 = F(fid, 1);
			ev2 = F(fid, 2);
			draws_highlighted_edge = true;
		}
		if (bc(1) < thresh)
		{
			ev1 = F(fid, 0);
			ev2 = F(fid, 2);
			draws_highlighted_edge = true;
		}
		if (bc(2) < thresh)
		{
			ev1 = F(fid, 0);
			ev2 = F(fid, 1);
			draws_highlighted_edge = true;
		}
		if (draws_highlighted_edge)
		{
			if (ev1 > ev2)
			{ // force ev1 < ev2, just for easier edge lookup for highlight
				swap(ev1, ev2);
			}

			if (prev_ev.first == ev1 && prev_ev.second == ev2)
				return true;

			// remove old hover
			if (prev_ev.first != -1)
			{ // there was an active edge, which is not active anymore, so delete it
				remove(highlighted_3d_edges, prev_ev);
				remove(highlighted_2d_edges, prev_uv1);
				if (prev_uv2.first != -1)
					remove(highlighted_2d_edges, prev_uv2);
			}

			num_uv_edges_found = find_corresponding_uv_edges(ev1, ev2);
			
			highlighted_3d_edges.push_back(pair<int, int>(ev1, ev2));
			highlighted_2d_edges.push_back(pair<int, int>(uv_edges[0].first, uv_edges[0].second));

			prev_ev = { ev1, ev2 };
			prev_uv1 = { uv_edges[0].first, uv_edges[0].second };

			if (uv_edges.size() == 2)
			{
				highlighted_2d_edges.push_back(pair<int, int>(uv_edges[1].first, uv_edges[1].second));
				prev_uv2 = { uv_edges[1].first, uv_edges[1].second };
			}
			update_colors = true;
			return true;
		}
	}
	else
	{
		mouse_over_3d_mesh = false;
	}

	if (prev_ev.first != -1)
	{ // there was an active edge, which is not active anymore, so delete it
		remove(highlighted_3d_edges, prev_ev);
		remove(highlighted_2d_edges, prev_uv1);
		if (prev_uv2.first != -1)
			remove(highlighted_2d_edges, prev_uv2);

		prev_ev = { -1, -1 };
		prev_uv1 = { -1, -1 };
		prev_uv2 = { -1, -1 };

		update_colors = true;
	}

	return false;
}

bool SolverPlugin::close_to_2d_mesh_edge()
{
	int fid;
	Eigen::Vector3f bc;
	viewer->core.viewport << 0, 0, 1200, 1350;
	double x = viewer->current_mouse_x;
	double y = viewer->core.viewport(3) - viewer->current_mouse_y;
	double thresh = 0.25;
	draws_highlighted_edge = false;
	int uv_ev1, uv_ev2;
	
	viewer->get_mesh(uv_id).dirty |= viewer->get_mesh(uv_id).DIRTY_OVERLAY_LINES;
	if (igl::unproject_onto_mesh(Vec2f(x, y), viewer->core.view * viewer->get_mesh(uv_id).model,
		viewer->core.proj, viewer->core.viewport, viewer->get_mesh(uv_id).V, viewer->get_mesh(uv_id).F, fid, bc))
	{
		mouse_over_2d_mesh = true;
		MatX3 viewer_uv_V = viewer->get_mesh(uv_id).V;
		if (bc(0) < thresh)
		{
			uv_ev1 = solver_wrapper->solver->Fs(fid, 1);
			uv_ev2 = solver_wrapper->solver->Fs(fid, 2);
			draws_highlighted_edge = true;
		}
		if (bc(1) < thresh)
		{
			uv_ev1 = solver_wrapper->solver->Fs(fid, 0);
			uv_ev2 = solver_wrapper->solver->Fs(fid, 2);
			draws_highlighted_edge = true;
		}
		if (bc(2) < thresh)
		{
			uv_ev1 = solver_wrapper->solver->Fs(fid, 0);
			uv_ev2 = solver_wrapper->solver->Fs(fid, 1);
			draws_highlighted_edge = true;
		}
		if (draws_highlighted_edge)
		{
			find_corresponding_mesh_edge(uv_ev1, uv_ev2, ev1, ev2);

			if (ev1 > ev2)
			{ // force ev1 < ev2, just for easier edge lookup for highlight
				swap(ev1, ev2);
			}

			if (prev_ev.first == ev1 && prev_ev.second == ev2)
				return true;

			// remove old hover
			if (prev_ev.first != -1)
			{ // there was an active edge, which is not active anymore, so delete it
				remove(highlighted_3d_edges, prev_ev);
				remove(highlighted_2d_edges, prev_uv1);
				if (prev_uv2.first != -1)
					remove(highlighted_2d_edges, prev_uv2);
			}

			num_uv_edges_found = find_corresponding_uv_edges(ev1, ev2);

			highlighted_3d_edges.push_back(pair<int, int>(ev1, ev2));
			highlighted_2d_edges.push_back(pair<int, int>(uv_edges[0].first, uv_edges[0].second));

			prev_ev = { ev1, ev2 };
			prev_uv1 = { uv_edges[0].first, uv_edges[0].second };

			if (uv_edges.size() == 2)
			{
				highlighted_2d_edges.push_back(pair<int, int>(uv_edges[1].first, uv_edges[1].second));
				prev_uv2 = { uv_edges[1].first, uv_edges[1].second };
			}
			update_colors = true;
			return true;
		}
	}
	else
	{
		mouse_over_2d_mesh = false;
	}

	if (prev_ev.first != -1)
	{ // there was an active edge, which is not active anymore, so delete it
		remove(highlighted_3d_edges, prev_ev);
		remove(highlighted_2d_edges, prev_uv1);
		if (prev_uv2.first != -1)
			remove(highlighted_2d_edges, prev_uv2);

		prev_ev = { -1, -1 };
		prev_uv1 = { -1, -1 };
		prev_uv2 = { -1, -1 };
	}

	return false;
}

bool SolverPlugin::close_to_2d_face()
{
	int fid;
	Eigen::Vector3f bc;
	viewer->core.viewport << 0, 0, 1200, 1350;
	double x = viewer->current_mouse_x;
	double y = viewer->core.viewport(3) - viewer->current_mouse_y;
	viewer->get_mesh(uv_id).dirty |= viewer->get_mesh(uv_id).DIRTY_OVERLAY_LINES;
	if (igl::unproject_onto_mesh(Vec2f(x, y), viewer->core.view * viewer->get_mesh(uv_id).model,
		viewer->core.proj, viewer->core.viewport, viewer->get_mesh(uv_id).V, viewer->get_mesh(uv_id).F, fid, bc))
	{
		mouse_over_2d_mesh = true;
		hovered_triangle = fid;
		update_colors = true;
		last_iteration_was_on_mesh = true;
		return true;
	}
	else
	{
		mouse_over_2d_mesh = false;
		if (last_iteration_was_on_mesh)
			update_colors = true;
		last_iteration_was_on_mesh = false;
	}
	hovered_triangle = -1;
	return false;
}

bool SolverPlugin::close_to_3d_face()
{
	int fid;
	Eigen::Vector3f bc;
	viewer->core.viewport << 1200, 0, 1200, 1350; // not needed with highlight viewport
	double x = viewer->current_mouse_x;
	double y = viewer->core.viewport(3) - viewer->current_mouse_y;
	viewer->get_mesh(mesh_id).dirty |= viewer->get_mesh(mesh_id).DIRTY_OVERLAY_LINES;
	if (igl::unproject_onto_mesh(Vec2f(x, y), viewer->core.view * viewer->get_mesh(mesh_id).model,
		viewer->core.proj, viewer->core.highlight_viewport, viewer->get_mesh(mesh_id).V, viewer->get_mesh(mesh_id).F, fid, bc))
	{
		mouse_over_3d_mesh = true;
		hovered_triangle = fid;
		update_colors = true;
		last_iteration_was_on_mesh = true;
		return true;
	}
	else
	{
		mouse_over_3d_mesh = false;
		if (last_iteration_was_on_mesh)
			update_colors = true;
		last_iteration_was_on_mesh = false;
	}
	hovered_triangle = -1;
	return false;
}

bool SolverPlugin::close_to_3d_vertex()
{
	int fid;
	Eigen::Vector3f bc;
	viewer->core.viewport << 1200, 0, 1200, 1350; // not needed with highlight viewport
	double x = viewer->current_mouse_x;
	double y = viewer->core.viewport(3) - viewer->current_mouse_y;
	viewer->get_mesh(mesh_id).dirty |= viewer->get_mesh(mesh_id).DIRTY_OVERLAY_LINES;
	if (igl::unproject_onto_mesh(Vec2f(x, y), viewer->core.view * viewer->get_mesh(mesh_id).model,
		viewer->core.proj, viewer->core.highlight_viewport, viewer->get_mesh(mesh_id).V, viewer->get_mesh(mesh_id).F, fid, bc))
	{
		mouse_over_3d_mesh = true;
		int idx;
		bc.maxCoeff(&idx);
		hovered_vertex = viewer->get_mesh(mesh_id).F(fid, idx);
		update_colors = true;
		last_iteration_was_on_mesh = true;
		return true;
	}
	else
	{
		mouse_over_3d_mesh = false;
		if (last_iteration_was_on_mesh)
			update_colors = true;
		last_iteration_was_on_mesh = false;
	}
	hovered_vertex = -1;
	return false;
}

void SolverPlugin::draw_dot_on_mouse_location()
{
	update_dot_on_mouse_pos();
	update_colors = true;
}

void SolverPlugin::update_dot_on_mouse_pos()
{
	double x = viewer->current_mouse_x;
	double y = viewer->core.viewport(3) - viewer->current_mouse_y;
	viewer->core.viewport << 0, 0, 1200, 1350;
	igl::unproject(RVec3(x, y, 0.), (viewer->core.view * viewer->get_mesh(uv_id).model).eval(), viewer->core.proj, viewer->core.viewport, dot_on_mouse_pos);
}

bool SolverPlugin::check_if_intersects_with_a_box_corner(int& ibox, int& icorner)
{
	update_dot_on_mouse_pos();
	RVec2 proj_curr_mouse;
	proj_curr_mouse << dot_on_mouse_pos(0), dot_on_mouse_pos(1);
	for (int i = 0; i < bbox_corners.size(); ++i)
	{
		if ((proj_curr_mouse - bbox_corners[i].first).norm() < 0.1)
		{
			ibox = i;
			icorner = 1;
			return true;
		}
		else if ((proj_curr_mouse - bbox_corners[i].second).norm() < 0.1)
		{
			ibox = i;
			icorner = 2;
			return true;
		}
	}
	return false;
}

void SolverPlugin::add_box()
{
	RVec2 A = bbox_corners.back().first;
	RVec2 B = bbox_corners.back().second;

	double xmin = min(A(0), B(0));
	double xmax = max(A(0), B(0));
	double ymin = min(A(1), B(1));
	double ymax = max(A(1), B(1));

	// for drawing set these as well
	top_left.push_back(RVec3(xmin, ymax, 0.));
	bottom_left.push_back(RVec3(xmin, ymin, 0.));
	bottom_right.push_back(RVec3(xmax, ymin, 0.));
	top_right.push_back(RVec3(xmax, ymax, 0.));

	solver_wrapper->get_slot();
	MatX2 X;
	solver_wrapper->solver->get_mesh(X);
	solver_wrapper->solver->energy->bbox->add_box(X, xmin, xmax, ymin, ymax);
	solver_wrapper->release_slot();

	bbox_exists = true;
}

void SolverPlugin::update_box(int boxid)
{
	RVec2 A = bbox_corners[boxid].first;
	RVec2 B = bbox_corners[boxid].second;

	double xmin = min(A(0), B(0));
	double xmax = max(A(0), B(0));
	double ymin = min(A(1), B(1));
	double ymax = max(A(1), B(1));

	// for drawing set these as well
	top_left[boxid] = RVec3(xmin, ymax, 0.);
	bottom_left[boxid] = RVec3(xmin, ymin, 0.);
	bottom_right[boxid] = RVec3(xmax, ymin, 0.);
	top_right[boxid] = RVec3(xmax, ymax, 0.);

	solver_wrapper->get_slot();
	solver_wrapper->solver->energy->bbox->update_box(boxid, xmin, xmax, ymin, ymax);
	solver_wrapper->release_slot();
}

void SolverPlugin::store_result(string fname)
{
	using namespace svg;

	string filename;
	if (fname.compare("") == 0)
		filename = igl::file_dialog_save();
	else
		filename = fname;

	// Store Scene
	save_scene(filename + ".scene");

	// from 0.5 -> 3
	viewer->get_mesh(uv_id).line_width *= 6.;

	// create SVG
	MatX3 VV = viewer->get_mesh(uv_id).V;
	MatX2 X = VV.block(0, 0, VV.rows(), 2);
	MatX3i XF = viewer->get_mesh(uv_id).F;

	// place bottom left vertex to origina of picture
	X.col(0).array() -= X.col(0).minCoeff();
	X.col(1).array() -= X.col(1).minCoeff();

	// set mean of svg vertices to 100
	double factor = 100. / X.mean();

	// find dimensions of whole mesh
	double xmin = factor * X.col(0).minCoeff();
	double xmax = factor * X.col(0).maxCoeff();
	double ymin = factor * X.col(1).minCoeff();
	double ymax = factor * X.col(1).maxCoeff();

	Dimensions dim(xmax - xmin, ymax - ymin);
	Document doc(filename + ".svg", svg::Layout(dim, svg::Layout::BottomLeft));

	RVec3i verts;
	Mat32 face;
	RVec2 v;
	for (int i = 0; i < XF.rows(); ++i)
	{
		verts = XF.row(i);
		igl::slice(X, verts, 1, face);
		doc << (svg::Polygon(svg::Color(svg::Color::Defaults::Transparent), Stroke(.5, svg::Color(0, 0, 0)))
			<< Point(factor * face(0, 0), factor * face(0, 1))
			<< Point(factor * face(1, 0), factor * face(1, 1))
			<< Point(factor * face(2, 0), factor * face(2, 1)));
	}
	doc.save();

	// create PNG
	// Allocate temporary buffers for 1200x1350 image
	MatC R(1200, 1350);
	MatC G(1200, 1350);
	MatC B(1200, 1350);
	MatC A(1200, 1350);

	viewer->core.camera_zoom = viewer->core.uv_camera_zoom;

	// Draw the scene in the buffers
	viewer->core.draw_buffer(viewer->get_mesh(uv_id), viewer->opengl[uv_id], true, R, G, B, A);

	// cutoff the 0 columns and rows to perfectly match the svg
	int start_row, end_row, start_col, end_col;
	find_mesh_boundaries_in_png_data(R, start_row, end_row, start_col, end_col);

	start_col *= resolution_factor;
	end_col *= resolution_factor;
	start_row *= resolution_factor;
	end_row *= resolution_factor;

	int w = end_col - start_col;
	int h = end_row - start_row;

	MatC Rbig(1200 * resolution_factor, 1350 * resolution_factor);
	MatC Gbig(1200 * resolution_factor, 1350 * resolution_factor);
	MatC Bbig(1200 * resolution_factor, 1350 * resolution_factor);
	MatC Abig(1200 * resolution_factor, 1350 * resolution_factor);

	// Draw the scene in the buffers
	viewer->core.draw_buffer(viewer->get_mesh(uv_id), viewer->opengl[uv_id], true, Rbig, Gbig, Bbig, Abig);

	// write 'with-boarder' png
	MatC R2 = Rbig.block(start_row, start_col, h, w);
	MatC G2 = Gbig.block(start_row, start_col, h, w);
	MatC B2 = Bbig.block(start_row, start_col, h, w);
	MatC A2 = Abig.block(start_row, start_col, h, w);

	igl::png::writePNG(R2, G2, B2, A2, filename + "_wb.png");

	// deactivate jaggy viewer lines
	viewer->get_mesh(uv_id).show_lines = false;
	// Draw the scene in the buffers
	viewer->core.draw_buffer(viewer->get_mesh(uv_id), viewer->opengl[uv_id], true, Rbig, Gbig, Bbig, Abig);

	R2 = Rbig.block(start_row, start_col, h, w);
	G2 = Gbig.block(start_row, start_col, h, w);
	B2 = Bbig.block(start_row, start_col, h, w);
	A2 = Abig.block(start_row, start_col, h, w);

	// Save it to a PNG
	igl::png::writePNG(R2, G2, B2, A2, filename + ".png");
	// reactivate jaggy viewer lines
	viewer->get_mesh(uv_id).show_lines = true;
	viewer->get_mesh(uv_id).line_width /= 6.;

	if (store_3d_mesh)
	{
		float old_lighting = viewer->core.lighting_factor;
		viewer->core.lighting_factor = 1.0f;
		viewer->core.camera_zoom = viewer->core.mesh_camera_zoom;

		// from 0.5 -> 3
		viewer->get_mesh(mesh_id).line_width *= 6.;

		// create PNG
		// Draw the scene in the buffers
		viewer->core.draw_buffer(viewer->get_mesh(mesh_id), viewer->opengl[mesh_id], true, R, G, B, A);

		// cutoff the 0 columns and rows to perfectly match the svg
		find_mesh_boundaries_in_png_data(R, start_row, end_row, start_col, end_col);

		start_col *= resolution_factor;
		end_col *= resolution_factor;
		start_row *= resolution_factor;
		end_row *= resolution_factor;

		w = end_col - start_col;
		h = end_row - start_row;

		// Draw the scene in the buffers
		viewer->core.draw_buffer(viewer->get_mesh(mesh_id), viewer->opengl[mesh_id], true, Rbig, Gbig, Bbig, Abig);

		// write 'with-boarder' png
		R2 = Rbig.block(start_row, start_col, h, w);
		G2 = Gbig.block(start_row, start_col, h, w);
		B2 = Bbig.block(start_row, start_col, h, w);
		A2 = Abig.block(start_row, start_col, h, w);

		igl::png::writePNG(R2, G2, B2, A2, filename + "_3d_wb.png");

		viewer->get_mesh(mesh_id).line_width /= 6.;
		viewer->core.lighting_factor = old_lighting;
	}
}

void SolverPlugin::find_mesh_boundaries_in_png_data(const MatC& R, int& start_row, int& end_row, int& start_col, int& end_col)
{
	// find start_col
	int i = 0;
	while (true)
	{
		int tot = R.col(i).sum();
		//if (R.col(i).sum() != 80 * resolution_factor)
		if (R.col(i).sum() != 80)
		{
			start_col = i;
			break;
		}
		++i;
	}
	// find end_col
	i = 1349;
	while (true)
	{
		if (R.col(i).sum() != 80)
		{
			end_col = i;
			break;
		}
		--i;
	}
	// find start_row
	i = 0;
	while (true)
	{
		if (R.row(i).sum() != 186)
		{
			start_row = i;
			break;
		}
		++i;
	}
	// find end_row
	i = 1199;
	while (true)
	{
		if (R.row(i).sum() != 186)
		{
			end_row = i;
			break;
		}
		--i;
	}
}

void SolverPlugin::export_mesh_to_svg_colored()
{
	ofstream out;
	out.open("test.svg", ofstream::out);
	if (!out.is_open())
	{
		cout << "cant open output file!" << endl;
		return;
	}

	MatX2 X;
	solver_wrapper->solver->get_mesh(X);

	MatX3i XF = solver_wrapper->solver->Fs;

	// place bottom left vertex to origina of picture
	X.col(0).array() -= X.col(0).minCoeff();
	X.col(1).array() -= X.col(1).minCoeff();

	double factor = 100.;

	// find dimensions of whole mesh
	double xmin = X.col(0).minCoeff();
	double xmax = X.col(0).maxCoeff();
	double ymin = X.col(1).minCoeff();
	double ymax = X.col(1).maxCoeff();

	// write header
	out << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.200000\" width=\"100\%\" height=\"100\%\" viewBox=\""
		<< xmin << " " << ymin << " " << xmax - xmin << " " << ymax - ymin << "\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">" << endl;
	
	for (int i = 0; i < XF.rows(); ++i)
	{
		RVec3i verts;
		Mat32 face;
		verts = XF.row(i);
		igl::slice(X, XF.row(i).eval(), 1, face);

		RVec2 v1 = face.row(0), v2 = face.row(1), v3 = face.row(2);

		RVec mid12 = (v1 + v2).array() / 2.;
		RVec mid23 = (v2 + v3).array() / 2.;
		RVec mid31 = (v3 + v1).array() / 2.;

		out << "\t<g>" << endl;
		out << "\t\t<defs>" << endl;

		// for each triangle we have 3 linear gradients, 2 paths and 1 filter in the defs
		out << "\t\t\t<linearGradient id=\"fadeA-" << i << "\" gradientUnits=\"userSpaceOnUse\" x1=\"" << v1(0) << "\" y1=\"" << v1(1) << "\" x2=\"" << mid23(0) << "\" y2=\"" << mid23(1) << "\">" << endl;
		out << "\t\t\t\t<stop offset = \"0\%\" stop-color=\"#FF0000\" />" << endl;
		out << "\t\t\t\t<stop offset = \"100\%\" stop-color=\"\#000000\" />" << endl;
		out << "\t\t\t</linearGradient>" << endl;

		out << "\t\t\t<linearGradient id=\"fadeB-" << i << "\" gradientUnits=\"userSpaceOnUse\" x1=\"" << v2(0) << "\" y1=\"" << v2(1) << "\" x2=\"" << mid31(0) << "\" y2=\"" << mid31(1) << "\">" << endl;
		out << "\t\t\t\t<stop offset = \"0\%\" stop-color=\"#00FF00\" />" << endl;
		out << "\t\t\t\t<stop offset = \"100\%\" stop-color=\"\#000000\" />" << endl;
		out << "\t\t\t</linearGradient>" << endl;

		out << "\t\t\t<linearGradient id=\"fadeC-" << i << "\" gradientUnits=\"userSpaceOnUse\" x1=\"" << v3(0) << "\" y1=\"" << v3(1) << "\" x2=\"" << mid12(0) << "\" y2=\"" << mid12(1) << "\">" << endl;
		out << "\t\t\t\t<stop offset = \"0\%\" stop-color=\"#0000FF\" />" << endl;
		out << "\t\t\t\t<stop offset = \"100\%\" stop-color=\"\#000000\" />" << endl;
		out << "\t\t\t</linearGradient>" << endl;

		out << "\t\t\t<path id=\"pathA-" << i << "\" d=\"M " << v1(0) << " " << v1(1) << " L " << v2(0) << " " << v2(1) << " " << v3(0) << " " << v3(1) << " Z\" fill=\"url(\#fadeA-" << i << ")\"/>" << endl;
		out << "\t\t\t<path id=\"pathB-" << i << "\" d=\"M " << v1(0) << " " << v1(1) << " L " << v2(0) << " " << v2(1) << " " << v3(0) << " " << v3(1) << " Z\" fill=\"url(\#fadeB-" << i << ")\"/>" << endl;

		out << "\t\t\t<filter id=\"Default" << i << "\">" << endl;
		out << "\t\t\t\t<feImage xlink:href=\"\#pathA-" << i << "\" result=\"layerA" << i << "\" x=\"0\" y=\"0\"/>" << endl;
		out << "\t\t\t\t<feImage xlink:href=\"\#pathB-" << i << "\" result=\"layerB" << i << "\" x=\"0\" y=\"0\"/>" << endl;
		out << "\t\t\t\t<feComposite in =\"layerA" << i << "\" in2=\"layerB" << i << "\" operator=\"arithmetic\" k1=\"0\" k2=\"1.0\" k3=\"1.0\" k4=\"0\" result=\"temp" << i << "\"/>" << endl;
		out << "\t\t\t\t<feComposite in =\"temp" << i << "\" in2=\"SourceGraphic\" operator=\"arithmetic\" k1=\"0\" k2=\"1.0\" k3=\"1.0\" k4=\"0\"/>" << endl;
		out << "\t\t\t</filter>" << endl;
		out << "\t\t</defs>" << endl;

		out << "\t\t<g stroke=\"none\" stroke-width=\"0\" shape-rendering=\"crispEdges\">" << endl;
		out << "\t\t\t<path d=\"M " << v1(0) << " " << v1(1) << " L " << v2(0) << " " << v2(1) << " " << v3(0) << " " << v3(1) << " Z\" fill=\"url(\#fadeC-" << i << ")\" filter=\"url(\#Default" << i << ")\"/>" << endl;
		out << "\t\t</g>" << endl;
		out << "\t</g>" << endl;
	}
	out << "</svg>" << endl;
	out.close();
}

void SolverPlugin::sample_down(int fac, MatC &M)
{
	int rows_orig = M.rows();
	int cols_orig = M.cols();

	int rows_new = rows_orig / 2;
	int cols_new = cols_orig / 2;

	// use fac*fac pixels to compute a single resulting pixel
	MatC Mt(rows_new, cols_new);

#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int col = 0; col < 2 * cols_new; col = col + 2)
	{
		Mat2C b;
		for (int row = 0; row < 2 * rows_new; row = row + 2)
		{
			b = M.block<2, 2>(row, col);
			int v1 = b(0, 0), v2 = b(0, 1), v3 = b(1, 0), v4 = b(1, 1);
			Mt(row / 2, col / 2) = (unsigned char)((v1 + v2 + v3 + v4) / 4);
		}
	}
	if (2 * cols_new < cols_orig)
	{
		VecC vc = M.col(cols_orig - 1);
		sample_down_vec(fac, vc);
		Mt.col(cols_new - 1) = vc;
	}
	if (2 * rows_new < rows_orig)
	{
		VecC vr = M.row(rows_orig - 1);
		sample_down_vec(fac, vr);
		Mt.row(rows_new - 1) = vr;
	}
	M = Mt;
}

void SolverPlugin::sample_down_vec(int fac, VecC &v)
{
	int rows_orig = v.rows();
	int rows_new = rows_orig / 2;

	VecC vt(rows_new);

	for (int i = 0; i < 2 * rows_new; i = i + 2)
	{
		VecC vb = v.block<2, 1>(i, 0);
		int v1 = vb(0), v2 = vb(1);
		vt(i / 2) = (unsigned char)((v1 + v2) / 2);
	}

	if (2 * rows_new < rows_orig)
		vt(rows_new - 1) = v(rows_orig - 1);

	v = vt;
}

void SolverPlugin::load_scene(const string& filename)
{
	ifstream in;
	in.open(filename, ifstream::in);
	if (!in.is_open())
	{
		cout << "cant open file" << endl;
		return;
	}

	int pos = filename.find_last_of('\\');
	mesh_filename = filename.substr(pos + 1, filename.length() - 6 - pos - 1) + ".obj";
	
	bool full_init = true;
	full_init &= read(solver_wrapper->solver->Fs, "Fs", in);
	full_init &= read(solver_wrapper->solver->Vs, "Vs", in);
	full_init &= read(F, "F", in);
	full_init &= read(V, "V", in);
	full_init &= read(solver_wrapper->solver->energy->separation->delta, "delta", in);
	full_init &= read(solver_wrapper->solver->energy->lambda, "lambda", in);
	full_init &= read(painting_colors, "painting_colors", in);
	full_init &= read(painting_weights, "painting_weights", in);
	full_init &= read(viewer->core.mesh_camera_zoom, "mesh_camera_zoom", in);
	full_init &= read(viewer->core.uv_camera_zoom, "uv_camera_zoom", in);
	full_init &= read(solver_wrapper->solver->energy->separation->V2V, "V2V", in);
	full_init &= read(solver_wrapper->solver->energy->separation->EVvar1, "EVvar1", in);
	full_init &= read(solver_wrapper->solver->energy->separation->EVvar2, "EVvar2", in);

	for (int i = 0; i < 3; ++i)
	{
		MatX3 currV;
		bool b = read(currV, "saved_state" + to_string(i), in);
		if (b)
		{
			double zoom;
			read(zoom, "state_camera_zoom" + to_string(i), in);
			MatX3 currN;
			read(currN, "state_normals" + to_string(i), in);
			saved_states[i] = currV;
			state_camera_zooms[i] = zoom;
			state_normals[i] = currN;
			cout << "camera state " << i << " loaded" << endl;
		}
	}

	in.close();
	
	solver_wrapper->solver->full_init_done = full_init;
	init_after_scene_load();
}

void SolverPlugin::save_scene(const string& filename)
{
	ofstream out;
	out.open(filename, ofstream::out);
	if (!out.is_open())
	{
		cout << "cant open file" << endl;
		return;
	}

	write(solver_wrapper->solver->Fs, "Fs", out);
	write(solver_wrapper->solver->Vs, "Vs", out);
	write(F, "F", out);
	write(V, "V", out);
	write(solver_wrapper->solver->energy->separation->delta, "delta", out);
	write(solver_wrapper->solver->energy->lambda, "lambda", out);
	write(painting_colors, "painting_colors", out);
	write(painting_weights, "painting_weights", out);
	write(viewer->core.mesh_camera_zoom, "mesh_camera_zoom", out);
	write(viewer->core.uv_camera_zoom, "uv_camera_zoom", out);
	write(solver_wrapper->solver->energy->separation->V2V, "V2V", out);
	write(solver_wrapper->solver->energy->separation->EVvar1, "EVvar1", out);
	write(solver_wrapper->solver->energy->separation->EVvar2, "EVvar2", out);
	
	for (int i = 0; i < 3; ++i)
	{
		if (saved_states[i].size() > 0)
		{
			write(saved_states[i], "saved_state" + to_string(i), out);
			write(state_camera_zooms[i], "state_camera_zoom" + to_string(i), out);
			write(state_normals[i], "state_normals" + to_string(i), out);
		}
	}

	out.close();
}

bool SolverPlugin::read(float& val, const string& name, ifstream& in)
{
	double d;
	bool ret = read(d, name, in);
	val = d;
	return ret;
}

bool SolverPlugin::read(double& val, const string& name, ifstream& in)
{
	string line, l, finline;
	bool last_was_space = false;

	while (getline(in, l))
	{	// read line and replace all multi spaces by a single comma	
		string currname;
		int i = 0;
		while (l[i] != ' ')
		{
			currname.push_back(l[i++]);
		}
		if (name.compare(currname) != 0)
			continue;

		i++;
		string sval;
		while (i < l.length())
			sval.push_back(l[i++]);

		val = std::atof(sval.c_str());
		return true;
	}
	return false;
}

template <typename Derived>
bool SolverPlugin::read(Eigen::PlainObjectBase<Derived>& M, const string& name, ifstream& in)
{
	string line, l, finline;
	bool last_was_space = false;

	while (getline(in, l))
	{	// read line and replace all multi spaces by a single comma	
		string currname;
		int i = 0;
		while (l[i] != ' ')
		{
			currname.push_back(l[i++]);
		}
		if (name.compare(currname) != 0)
			continue;

		line = l;
		int cnt = 0;
		for (int i = 0; i < l.length(); ++i)
		{
			if (l[i] == ' ')
			{
				if (!last_was_space)
				{
					line[cnt++] = ',';
				}
				last_was_space = true;
			}
			else
			{
				line[cnt++] = l[i];
				last_was_space = false;
			}
		}
		finline = line.substr(0, cnt);
		finline.push_back(',');
		// read the values
		vector<double> vals;
		double val;
		string srows;
		int iidx = currname.length() + 1;
		while (finline[iidx] != ',')
		{
			srows.push_back(finline[iidx++]);
		}
		int rows = std::atof(srows.c_str());
		iidx++;
		string scols;
		while (finline[iidx] != ',')
		{
			scols.push_back(finline[iidx++]);
		}
		int cols = std::atof(scols.c_str());
		iidx++;
		int ii = iidx;
		while (ii < finline.length())
		{
			string snum;
			while (finline[ii] != ',')
				snum.push_back(finline[ii++]);
			ii++;
			val = std::atof(snum.c_str());
			vals.push_back(val);
		}
		M = Eigen::PlainObjectBase<Derived>::Zero(rows, cols);
		int vptr = 0;
		for (int c = 0; c < cols; ++c)
		{
			for (int r = 0; r < rows; ++r)
			{
				M(r, c) = vals[vptr++];
			}
		}
		return true;
	}
	return false;
}

bool SolverPlugin::read(SpMat& M, const string& name, ifstream& in)
{
	string line, l, finline;
	bool last_was_space = false;

	while (getline(in, l))
	{	// read line and replace all multi spaces by a single comma	
		string currname;
		int i = 0;
		while (l[i] != ' ')
		{
			currname.push_back(l[i++]);
		}
		if (name.compare(currname) != 0)
			continue;

		line = l;
		int cnt = 0;
		for (int i = 0; i < l.length(); ++i)
		{
			if (l[i] == ' ')
			{
				if (!last_was_space)
				{
					line[cnt++] = ',';
				}
				last_was_space = true;
			}
			else
			{
				line[cnt++] = l[i];
				last_was_space = false;
			}
		}
		finline = line.substr(0, cnt);
		// read the values
		vector<double> vals;
		double val;
		string srows;
		int iidx = currname.length() + 1;
		while (finline[iidx] != ',')
		{
			srows.push_back(finline[iidx++]);
		}
		int rows = std::atof(srows.c_str());
		iidx++;
		string scols;
		while (finline[iidx] != ',')
		{
			scols.push_back(finline[iidx++]);
		}
		int cols = std::atof(scols.c_str());
		iidx++;
		int ii = iidx;
		while (ii < finline.length())
		{
			string snum;
			while (finline[ii] != ',')
				snum.push_back(finline[ii++]);
			ii++;
			val = std::atof(snum.c_str());
			vals.push_back(val);
		}
		M = SpMat(rows, cols);
		Triplets triplets;
		for (int c = 0; c < vals.size(); c += 3)
		{
			triplets.push_back(Tripletd(vals[c], vals[c + 1], vals[c + 2]));
		}
		M.setFromTriplets(triplets.begin(), triplets.end());
		return true;
	}
	return false;
}

void SolverPlugin::find_largest_non_overlapping_set()
{
	nonOverlappingTriIDs.clear();
	patchBoundaryEdges.clear();
	// find clicked triangle
	int fid;
	Eigen::Vector3f bc;
	viewer->core.viewport << 0, 0, 1200, 1350;
	double x = viewer->current_mouse_x;
	double y = viewer->core.viewport(3) - viewer->current_mouse_y;
	MatX3 solver_Vs_zero_extended = MatX3::Zero(solver_wrapper->solver->Vs.rows(), 3);
	solver_Vs_zero_extended.block(0, 0, solver_wrapper->solver->Vs.rows(), 2) = solver_wrapper->solver->Vs;
	if (!igl::unproject_onto_mesh(Vec2f(x, y), viewer->core.view * viewer->get_mesh(uv_id).model,
		viewer->core.proj, viewer->core.viewport, solver_Vs_zero_extended, solver_wrapper->solver->Fs, fid, bc))
		return;
	
	MatX3 snap_UV;
	MatX3i Fsnap;
	// pad uv with zeros
	MatX2 UV;

	// Snap UVs
	solver_wrapper->solver->get_mesh(UV);
	MatX3 padded_UV = MatX3::Zero(UV.rows(), 3);
	padded_UV.block(0, 0, UV.rows(), 2) = UV;
	auto& Fs = solver_wrapper->solver->Fs;
	snap_soup_mat(padded_UV, Fs, snap_UV, Fsnap);
	map<pair<int, int>, pair<int, int>> snapEdge2soupEdge;

	auto MakeSortedPair = [](int i, int j) { return (i < j) ? make_pair(i, j) : make_pair(j, i); };

	for (int i = 0; i < Fs.rows(); i++)
	{
		snapEdge2soupEdge[MakeSortedPair(Fsnap(i, 0), Fsnap(i, 1))] = make_pair(Fs(i, 0), Fs(i, 1));
		snapEdge2soupEdge[MakeSortedPair(Fsnap(i, 1), Fsnap(i, 2))] = make_pair(Fs(i, 1), Fs(i, 2));
		snapEdge2soupEdge[MakeSortedPair(Fsnap(i, 2), Fsnap(i, 0))] = make_pair(Fs(i, 2), Fs(i, 0));
	}

//  Based on Roi 2017
		// Get largest nonoverlappig set
		queue<int> nextTriangles;
		vector<bool> visitedTriangles(Fsnap.rows(), false);
		set<pair<int, int>> snapPatchBoundaryEdges;
		
		auto AddOrRemoveEdge = [&snapPatchBoundaryEdges](pair<int, int> p) {
			if (snapPatchBoundaryEdges.count(p) == 0)
				snapPatchBoundaryEdges.insert(p);
			else
				snapPatchBoundaryEdges.erase(p);
		};
	
		nextTriangles.push(fid);
		snapPatchBoundaryEdges.insert(MakeSortedPair(Fsnap(fid, 0), Fsnap(fid, 1)));
		snapPatchBoundaryEdges.insert(MakeSortedPair(Fsnap(fid, 1), Fsnap(fid, 2)));
		snapPatchBoundaryEdges.insert(MakeSortedPair(Fsnap(fid, 2), Fsnap(fid, 0)));
	
		Eigen::MatrixX3i TT, TTi;
		igl::triangle_triangle_adjacency(Fsnap, TT, TTi);
		while (nextTriangles.size() > 0)
		{
			int id = nextTriangles.front();
			nextTriangles.pop();
			visitedTriangles[id] = true;
	
			nonOverlappingTriIDs.push_back(id);
			
			for(int i=0; i<3; i++)
			{
				int ntId = TT(id, i);
				if (ntId != -1 && visitedTriangles[ntId] == false)
				{
					visitedTriangles[ntId] = true;
					vector < pair<int, int>> newTriangle;
					newTriangle.push_back(MakeSortedPair(Fsnap(ntId, 0), Fsnap(ntId, 1)));
					newTriangle.push_back(MakeSortedPair(Fsnap(ntId, 1), Fsnap(ntId, 2)));
					newTriangle.push_back(MakeSortedPair(Fsnap(ntId, 2), Fsnap(ntId, 0)));
					//check intersections
					bool flag = false;
					for (auto itr1 = snapPatchBoundaryEdges.begin(); itr1 != snapPatchBoundaryEdges.end(); itr1++) {
						for (auto itr2 = newTriangle.begin(); itr2 != newTriangle.end(); itr2++)
						{
							auto e1(*itr1) , e2(*itr2);
							if (e1.first==e2.first || e1.first==e2.second ||
								e1.second == e2.first || e1.second == e2.second)
								continue;
							Vec2 p1 = snap_UV.row(e1.first).leftCols(2);
							Vec2 p2 = snap_UV.row(e1.second).leftCols(2);
							Vec2 q1 = snap_UV.row(e2.first).leftCols(2);
							Vec2 q2 = snap_UV.row(e2.second).leftCols(2);
							double t, u;

							if(testIntersection(p1,p2,q1,q2))
								flag = true;
							if (flag)
								break;
						}
						if (flag)
							break;
					}
// 					edges ee;
// 					for (auto edge : boundaryEdges)
// 					{
// 						point p1 = { snap_UV(edge.first,0), snap_UV(edge.first, 1) };
// 						point p2 = { snap_UV(edge.second,0), snap_UV(edge.second, 1) };
// 						ee.insert(make_pair(p1, p2));
// 					}
// 					bool flag = is_simple(ee);
					if (!flag)
					{
						nextTriangles.push(ntId);
						AddOrRemoveEdge(MakeSortedPair(Fsnap(ntId, 0), Fsnap(ntId, 1)));
						AddOrRemoveEdge(MakeSortedPair(Fsnap(ntId, 1), Fsnap(ntId, 2)));
						AddOrRemoveEdge(MakeSortedPair(Fsnap(ntId, 2), Fsnap(ntId, 0)));
					}
				}
			}
		}
		Eigen::MatrixX2i E(snapPatchBoundaryEdges.size(), 2);
		Eigen::VectorXi B(snapPatchBoundaryEdges.size());
		int ind = 0;
		for (auto e : snapPatchBoundaryEdges)
		{
			E(ind, 0) = e.first;
			E(ind++, 1) = e.second;
		}
// 		igl::is_boundary_edge(E, Fsnap, B);
// 		for (int i = 0; i < B.size(); i++)
// 			if (B(i))
// 				snapPatchBoundaryEdges.erase(make_pair(E(i, 0), E(i, 1)));

		for (auto e : snapPatchBoundaryEdges)
			patchBoundaryEdges.insert(snapEdge2soupEdge[e]);
}

void SolverPlugin::freePatchBoundaryEdges()
{
	int uv_ev1, uv_ev2;
	for (auto e: patchBoundaryEdges)
	{
		// uv_ev is one edge
		uv_ev1 = e.first;
		uv_ev2 = e.second;
		// ev are on the original mesh
		find_corresponding_mesh_edge(uv_ev1, uv_ev2, ev1, ev2);

		num_uv_edges_found = find_corresponding_uv_edges(ev1, ev2);

		react_to_edge_click();
	}
}

bool SolverPlugin::testIntersection(Vec2 p1, Vec2 p2, Vec2 q1, Vec2 q2)
{
	auto cross2d = [](Vec2 a0, Vec2 a1) {
		return a0(0)*a1(1) - a0(1)*a1(0);
	};

	auto differentSides = [&](Vec2 a0, Vec2 a1, Vec2 o, Vec2 d) {
		return (cross2d(d, (a0 - o).eval()) * cross2d(d, (a1 - o).eval()) <= 0);
	};

	return differentSides(p1, p2, q1, q2 - q1) && differentSides(q1, q2, p1, p2 - p1);
}

inline void SolverPlugin::write(const Mat& M, const string& name, ofstream& out)
{
	RVec v;
	to_rvec(M, v);
	out << name << " " << M.rows() << " " << M.cols() << " " << v << endl;
}

inline void SolverPlugin::write(const MatX3i& M, const string& name, ofstream& out)
{
	RVeci v;
	to_rvec(M, v);
	out << name << " " << M.rows() << " " << M.cols() << " " << v << endl;
}

inline void SolverPlugin::write(const MatX2& M, const string& name, ofstream& out)
{
	RVec v;
	to_rvec(M, v);
	out << name << " " << M.rows() << " " << M.cols() << " " << v << endl;
}

inline void SolverPlugin::write(const MatX3& M, const string& name, ofstream& out)
{
	RVec v;
	to_rvec(M, v);
	out << name << " " << M.rows() << " " << M.cols() << " " << v << endl;
}

inline void SolverPlugin::write(const SpMat& M, const string& name, ofstream& out)
{
	out << name << " " << M.rows() << " " << M.cols() << " ";
	for (int i = 0; i < M.outerSize(); ++i)
	{
		for (SpMat::InnerIterator it(M, i); it; ++it)
		{
			out << it.row() << " " << it.col() << " " << it.value() << " ";
		}
	}
	out << endl;
}

inline void SolverPlugin::write(double v, const string& name, ofstream& out)
{
	out << name << " " << v << endl;
}

inline void SolverPlugin::to_rvec(const Mat& M, RVec& v)
{
	v = Eigen::Map<const RVec>(M.data(), 1, M.rows()*M.cols());
}

inline void SolverPlugin::to_rvec(const MatX3i& M, RVeci& v)
{
	v = Eigen::Map<const RVeci>(M.data(), 1, M.rows()*M.cols());
}

inline void SolverPlugin::to_rvec(const MatX2& M, RVec& v)
{
	v = Eigen::Map<const RVec>(M.data(), 1, M.rows()*M.cols());
}

inline void SolverPlugin::to_rvec(const MatX3& M, RVec& v)
{
	v = Eigen::Map<const RVec>(M.data(), 1, M.rows()*M.cols());
}

void SolverPlugin::set_hit_triangle_as_active()
{
	int fid;
	Eigen::Vector3f bc;
	viewer->core.viewport << 0, 0, 1200, 1350;
	double x = viewer->current_mouse_x;
	double y = viewer->core.viewport(3) - viewer->current_mouse_y;
	last_hit_triangle = hit_triangle;
	MatX3 solver_Vs_zero_extended = MatX3::Zero(solver_wrapper->solver->Vs.rows(), 3);
	solver_Vs_zero_extended.block(0, 0, solver_wrapper->solver->Vs.rows(), 2) = solver_wrapper->solver->Vs;
	if (igl::unproject_onto_mesh(Vec2f(x, y), viewer->core.view * viewer->get_mesh(uv_id).model,
		viewer->core.proj, viewer->core.viewport, solver_Vs_zero_extended, solver_wrapper->solver->Fs, fid, bc))
	{
		hit_triangle = fid;
		if (!solver_wrapper->if_its_a_fixed_triangle_remove_and_store(hit_triangle, triangle_pos_down))
			igl::slice(solver_wrapper->solver->Vs, solver_wrapper->solver->Fs.row(fid).eval(), 1, triangle_pos_down);
		update_triangle_position(triangle_pos_down);
		update_colors = true;
		move_triangle_enabled = true;
	}
	else
		hit_triangle = -1;
	solver_wrapper->set_active_triangle(hit_triangle);
}

void SolverPlugin::add_triangle_to_fixed_position_set()
{
	int fid;
	Eigen::Vector3f bc;
	viewer->core.viewport << 0, 0, 1200, 1350;
	double x = viewer->current_mouse_x;
	double y = viewer->core.viewport(3) - viewer->current_mouse_y;
	MatX3 solver_Vs_zero_extended = MatX3::Zero(solver_wrapper->solver->Vs.rows(), 3);
	solver_Vs_zero_extended.block(0, 0, solver_wrapper->solver->Vs.rows(), 2) = solver_wrapper->solver->Vs;
	if (igl::unproject_onto_mesh(Vec2f(x, y), viewer->core.view * viewer->get_mesh(uv_id).model,
		viewer->core.proj, viewer->core.viewport, solver_Vs_zero_extended, solver_wrapper->solver->Fs, fid, bc))
	{
		Mat32 pos;
		igl::slice(solver_wrapper->solver->Vs, solver_wrapper->solver->Fs.row(fid).eval(), 1, pos);
		solver_wrapper->add_or_remove_triangle_from_fixed_position_set(fid, pos);
		update_colors = true;
	}
}

void SolverPlugin::update_triangle_hit_color()
{
	if (last_hit_triangle == -1 && hit_triangle == -1)
		return;
	if (last_hit_triangle == hit_triangle)
	{
		move_triangle_enabled = true;
		return;
	}
	if (last_hit_triangle == -1)
	{
		move_triangle_enabled = true;
	}
	else if (hit_triangle == -1)
	{
		move_triangle_enabled = false;
	}
	else
	{
		move_triangle_enabled = true;
	}
}

void SolverPlugin::translate(double offset_x, double offset_y)
{
	MatX3 new_mesh_pos = mesh_pos_down;
	new_mesh_pos.col(0).array() += offset_x;
	new_mesh_pos.col(1).array() += offset_y;
	set_3d_soup_vertices_with_connected_vertices(new_mesh_pos);
	V = new_mesh_pos;
}

void SolverPlugin::translate_triangle(double offset_x, double offset_y)
{
	Mat32 new_pos = triangle_pos_down;
	RVec2 offset(offset_x, offset_y);
	new_pos.row(0) += offset;
	new_pos.row(1) += offset;
	new_pos.row(2) += offset;
	moving_triangle_pos = new_pos;
}

void SolverPlugin::translate_uv_mesh(double offset_x, double offset_y)
{
	MatX3 new_mesh_pos = uv_mesh_pos_down;
	new_mesh_pos.col(0).array() += offset_x;
	new_mesh_pos.col(1).array() += offset_y;
	viewer->get_mesh(uv_id).set_vertices(new_mesh_pos);
}

void SolverPlugin::rotate(double phi_x, double phi_y)
{
	Rx.block<2, 2>(1, 1) << cos(phi_x), sin(phi_x), -sin(phi_x), cos(phi_x);
	Ry.row(0) << cos(phi_y), 0, sin(phi_y);
	Ry.row(2) << -sin(phi_y), 0, cos(phi_y);
	Mat3 R = Ry*Rx;
	MatX3 Vtemp = mesh_pos_down;
	RVec3 col_mean = Vtemp.colwise().mean();
	Vtemp.rowwise() -= col_mean;
	MatX3 Vrot = Vtemp * R;
	Vrot.rowwise() += col_mean;
	set_3d_soup_vertices_with_connected_vertices(Vrot);
	viewer->get_mesh(0).set_normals(mesh_3d_normals_down * R);
	V = Vrot;
}

void SolverPlugin::scale_3d_mesh(bool up)
{
	MatX3 Vtemp = V;
	Vtemp.array() *= up ? 1.1 : 0.9;
	set_3d_soup_vertices_with_connected_vertices(Vtemp);
	V = Vtemp;
}

void SolverPlugin::init_rotation_matrices()
{
	Rx = Mat3::Zero();
	Rx(0, 0) = 1;
	Ry = Mat3::Zero();
	Ry(1, 1) = 1;
}

bool SolverPlugin::remove(map<pair<int, int>, MatX3>& l, const pair<int, int>& el)
{
	for (map<pair<int, int>, MatX3>::iterator it = l.begin(); it != l.end(); ++it)
	{
		if (it->first == el)
		{
			l.erase(it);
			return true;
		}
	}
	return false;
}

bool SolverPlugin::remove(list<pair<int, int>>& l, const pair<int, int>& el)
{
	for (auto it = l.begin(); it != l.end(); ++it)
	{
		if (*it == el)
		{
			l.erase(it);
			return true;
		}
	}
	return false;
}

void SolverPlugin::add_viewport_highlights(int x)
{ // total width = 2400
	if (!mesh_loaded || !update_highlights)
		return;

	viewer->get_mesh(vp_id).lines.resize(0, Eigen::NoChange);
	Mat P1 = Mat(4, 3);
	Mat P2 = Mat(4, 3);
	P1 << 2, 2.06, 0,
		-1.83, 2.05, 0,
		2, -2.06, 0,
		1.83, 2.05, 0;
	P2 << -2, 2.06, 0,
		-1.83, -2.1, 0,
		-2, -2.06, 0,
		1.83, -2.1, 0;
	if (x <= 1200)
	{ // uv/left side
		mouse_on_uv_side = true;
		viewer->core.highlight_viewport << 0, 0, 1200, 1350;
		viewer->core.highlight_camera_zoom = viewer->core.uv_camera_zoom;
	}
	else
	{ // 3d/right side
		mouse_on_uv_side = false;
		viewer->core.highlight_viewport << 1200, 0, 1200, 1350;
		viewer->core.highlight_camera_zoom = viewer->core.mesh_camera_zoom;
	}
	viewer->get_mesh(vp_id).add_edges(P1, P2, white);
}

void SolverPlugin::set_3d_soup_vertices_with_connected_vertices(const MatX3& newV)
{
#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int i = 0; i < F.rows(); ++i)
	{
		RVec3i verts;
		Mat3 face;
		verts = F.row(i);
		igl::slice(newV, verts, 1, face);
		mesh_soup_V.block<3, 3>(3 * i, 0) = face;
	}
	viewer->get_mesh(mesh_id).set_mesh(mesh_soup_V, mesh_soup_F);
}

void SolverPlugin::get_3d_mesh_connected_vertices(MatX3& Vout, MatX3i Fout)
{
	Vout = V;
	Fout = F;
}