#pragma once

#ifndef SOLVERPLUGIN_H
#define SOLVERPLUGIN_H

#include "ApplicationPipelineStep.h"
#include "DataObject.h"
#include "EigenTypes.h"
#include "SolverWrapper.h"

#include <thread>
#include <nanogui/slider.h>

#ifndef INT_INF
#define INT_INF numeric_limits<int>::max()
#endif

using namespace std;
using namespace nanogui;

class SolverPlugin : public ApplicationPipelineStep
{
public:
	enum class Param { LAMBDA, DELTA, BOUND, POSITION_WEIGHT };
	enum class Mode : int { MOVE = 0, EDGE_CUTTING, FACE_POSITIONING, PAINTING, VERTEX_CLICKING, BBOX_DRAWING };
	enum class PainterMode : int { FIX, MOVE, ERASE };

	SolverPlugin();

	void InitPipelineStep();
	bool Load(string filename);
	void add_slider();
	void add_color_clamp_slider(const string& name, const shared_ptr<double>& max_value, const shared_ptr<double>& value);
	void find_edge_lenghts_for_separation();
	void find_new_edges_after_triangulation();
	void snap_soup_mat(const MatX3& Ms, const MatX3i& Fs, MatX3& M, MatX3i& Md);
	void export_uv_to_obj();
	void read_quad_mesh(const string& quad_name);
	bool edge_exists_in_quad_mesh(int p1, int p2);
	int find_corresponding_pair_index(int i1, int i2);
	void add_texture_slider(nanogui::Window* window, double& var, const string& name);
	void init();
	void init_after_scene_load();
	
	bool MouseMove(int mouse_x, int mouse_y);
	void process_mouse_move();

	bool MouseDown(int button, int modifier);
	bool MouseUp(int button, int modifier);
	bool MouseScroll(float delta_y);

	void scale_3d_mesh(bool up);
	void scale_uv_mesh(bool up);
	void rotate(double phi_x, double phi_y);
	void translate(double offset_x, double offset_y);
	void translate_uv_mesh(double offset_x, double offset_y);
	void translate_triangle(double offset_x, double offset_y);
	
	bool PreDraw();
	bool KeyDown(int key, int modifiers);
	bool KeyUp(int key, int modifiers);

	void update_mesh();
	void start_solver_thread();
	void run_solver_and_store_every_iteration();
	void stop_solver_thread();
	void update_energy_param(Param p, double value, bool start_solver = true);
	void update_triangle_position(const Mat32& target_triangle_pos);

	void react_to_edge_click();
	void react_to_edge_click_force_merge();
	int find_corresponding_uv_edges(int v1, int v2);
	void find_corresponding_mesh_edge(int uv_ev1, int uv_ev2, int& ev1, int& ev2);

	void apply_weight_onto_hit_corners();

	void build_colored_squares();
	void update_boxes_colors();

	bool colored_squares_already_built = false;

	bool rotation_enabled = false;
	bool translation_enabled = false;
	bool move_triangle_enabled = false;
	bool uv_translation_enabled = false;
	bool applying_weight_enabled = false;

	int last_mouse_x, last_mouse_y;

	thread solver_thread;
	unique_ptr<SolverWrapper> solver_wrapper;

	// lambda slider
	Slider *slider;

	// vp = viewport
	// xh = crosshair
	// cs = colored squares
	unsigned int uv_id = 0, mesh_id = 0, vp_id = 0, xh_id = 0, cs_id = 0;

	Mode mode = Mode::FACE_POSITIONING;

	double max_weighting_val = 10000;
	double weighting_step = 1000;

	bool set_edge_lenghts_to_average = false;
	double texture_size = 0.5;
	double max_texture_val = 10.;

	Utils::Init uv_init = Utils::Init::ISOMETRIC;
	bool wants_to_use_harmonic_mapping = false;

	int solver_iter_counter = 0;
	string foldername = "harmonic";
	bool use_num_steps = false;

	int resolution_factor = 4;

	bool show_harmonic_cut = true;

	bool update_colors = false;

	bool store_3d_mesh = true;

	bool show_separation_error = true;
	bool show_distortion_error = false;

	shared_ptr<double> sep_color_clamp = make_shared<double>(0.5);
	shared_ptr<double> dist_color_clamp = make_shared<double>(0.5);

	shared_ptr<double> max_sep_color_value = make_shared<double>(3.0);
	shared_ptr<double> max_dist_color_value = make_shared<double>(1.0);

	string mesh_filename;

	double no_seam_force = 100.;

	bool save_state_active = false;
	vector<MatX3> saved_states = vector<MatX3>(3, MatX3());
	vector<double> state_camera_zooms = vector<double>(3);
	vector<MatX3> state_normals = vector<MatX3>(3, MatX3());

	bool show_fixed_3d_edge_highlights = true;
	bool show_quad_mesh = false;

	MatX3 Vq;
	MatX4i Fq;

	MatX3 old_normals;
	MatX3 face_N;

	bool load_uv_from_file = false;

	Vec3 xh_center_3d;

private:
	// The 3d mesh
	MatX3 V;
	MatX3i F;

	// the cut mesh
	MatX3 V_cut;
	MatX3i F_cut;

	// a representation of the 3d mesh as a soup
	// this is needed for painting on it, while being able
	// to translate each 2d soup triangle onto the
	// 3d soup and vice-versa.
	MatX3 mesh_soup_V;
	MatX3i mesh_soup_F;

	// mapping from original/connected 3d mesh, to the 3d soup mesh vertices
	// a vertex is split into either 2 or 3 soup vertices
	//map<int, vector<int>> mesh_map;
	map<int, vector<int>> mesh_map_to_soup, mesh_map_to_orig;

	// Rotation matrices
	Mat3 Rx, Ry;

	MatX3 uv_triangle_colors, mesh_triangle_colors;
	int hit_triangle = -1, last_hit_triangle = -1, hovered_triangle = -1, hovered_vertex = -1;
	//RVec3 last_hit_triangle_color;

	// when loading a mesh with cut option, place two points and map these to a circle (boundary)
	bool is_in_initial_cut_mode = false;

	// we need 2 points on the 3d mesh to place the boundary on as start and end of the cut
	pair<int, int> initial_cut = pair<int, int>(-1, -1);

	// path through vertices, initial cut
	vector<int> path;
	
	inline string removeTrailingZeros(string&& s);
	void init_rotation_matrices();
	void update_triangle_hit_color();
	void add_triangle_to_fixed_position_set();
	void set_hit_triangle_as_active();

	bool remove(map<pair<int, int>, MatX3>& l, const pair<int, int>& el);
	bool remove(list<pair<int, int>>& l, const pair<int, int>& el);

	bool close_to_3d_mesh_edge();
	bool close_to_2d_mesh_edge();

	bool close_to_2d_face();
	bool close_to_3d_face();
	bool close_to_3d_vertex();

	void reset_edge_cutting();
	void reset_face_positioning();
	void reset_painting();

	void add_viewport_highlights(int x);

	void show_menus();
	void hide_menus();

	void hide_all_menus();
	void show_all_menus();

	bool draws_highlighted_edge = false;
	// only force an update of the coloring, when we leave the mesh
	// and not when we move the mouse in empty space, there no
	// recoloring is needed (besides normal mesh updates from the solver ofc)
	bool last_iteration_was_on_mesh = false;

	void set_3d_soup_vertices_with_connected_vertices(const MatX3& newV);
	void get_3d_mesh_connected_vertices(MatX3& V, MatX3i F);
	void map_to_soup_matrix(const MatX2 in, const MatX3i& F, MatX2& out);

	// tria-tria-intersection stuff
	typedef pair<Vec3, Vec3> Edge;
	void check_intersections();
	bool line_intersects_triangle(const Edge& e, const Mat3& T);
	void find_edges_of_triangle(int face, const MatX3& Vs, const MatX3i& Fs, Edge& e1, Edge& e2, Edge& e3);
	bool test_intersection(const Edge& e1_1, const Edge& e2_1);
	void toggle_free_overlapping_faces();
	void toggle_free_face(const RVec3i& face);

	// bbox stuff
	void draw_dot_on_mouse_location();
	void update_dot_on_mouse_pos();
	bool check_if_intersects_with_a_box_corner(int& ibox, int& icorner);
	void add_box();
	void update_box(int boxid);

	void store_result(string filename = "");
	void find_mesh_boundaries_in_png_data(const MatC& R, int& start_row, int& end_row, int& start_col, int& end_col);
	void export_mesh_to_svg_colored();
	void sample_down(int fac, Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> &M);
	void sample_down_vec(int fac, Eigen::Matrix<unsigned char, Eigen::Dynamic, 1> &v);

	void load_scene(const string& filename);
	void save_scene(const string& filename);

	void write(const Mat& M, const string& name, ofstream& out);
	void write(const MatX3i& M, const string& name, ofstream& out);
	void write(const MatX2& M, const string& name, ofstream& out);
	void write(const MatX3& M, const string& name, ofstream& out);
	void write(const SpMat& M, const string& name, ofstream& out);
	void write(double v, const string& name, ofstream& out);

	void to_rvec(const Mat& M, RVec& v);
	void to_rvec(const MatX3i& M, RVeci& v);
	void to_rvec(const MatX2& M, RVec& v);
	void to_rvec(const MatX3& M, RVec& v);

	template <typename Derived>
	bool read(Eigen::PlainObjectBase<Derived>& M, const string& name, ifstream& in);
	bool read(SpMat& M, const string& name, ifstream& in);
	bool read(double& val, const string& name, ifstream& in);
	bool read(float& val, const string& name, ifstream& in);
	void run_sorkine2002();
	void Seperate_overlaps_using_sorkine2002();
	void find_largest_non_overlapping_set();
	void freePatchBoundaryEdges();
	bool testIntersection(Vec2 p1, Vec2 p2, Vec2 q1, Vec2 q2);
	bool show_intersections = false;
	bool free_overlapping_triangles = false;
	Vecb overlapping_triangles;

	Vec ones_vec;

	// Current position of the moving triangle
	Mat32 moving_triangle_pos;

	// the two vertices involved in the edge click (on the 3d mesh)
	int ev1, ev2;

	// in case of not hitting an edge, we have to remove the previously active (hovered) from the highlight set
	pair<int, int> prev_ev, prev_uv1, prev_uv2;

	// The UV vertex pair(s) mapped from the 3d mesh vertex pair ev1-ev2
	vector<pair<int, int>> uv_edges;

	// number of edges found when mapping from 3d to 2d
	// needed to clear highlighted edges (conservative resize)
	int num_uv_edges_found;

	// a map for all edges already clicked on (fixedly highlighted) to handle Esep
	// 3d edge -> 1 or 2 2d edges
	map<pair<int, int>, pair<int, int>> zeroed_esep_columns;
	map<pair<int, int>, pair<int, int>>::iterator zeroed_map_it;

	// a single lookup to see if an edge is already highlighted or not (edge = vi-vj, vi<vj)
	map<pair<int, int>, MatX3> highlighted_edges;

	// Last edge vertices of interest
	RVec3 P1, P2;
	RVec3 UVP1, UVP2, UVP3, UVP4;
	RVec3 C, C_hover, white, red, C_merge, zero3, ones3, black;

	// vertex-vertex pairs (== edges) which are currently highlighted in the 3d and uv/2d mesh
	list<pair<int, int>> highlighted_3d_edges, highlighted_2d_edges;

	// list of FIXED edges, which have been clicked on
	list<pair<int, int>> fixed_highlighted_3d_edges, fixed_highlighted_2d_edges, sep_highlights_3d;
	map<pair<int, int>, RVec3> fixed_highlight_3d_colors, fixed_highlight_2d_colors;

	Mat32 triangle_pos_down;

	MatX3 mesh_pos_down, uv_mesh_pos_down;
	RVec3 last_mouse_pos, projected_mouse_down, point_pos_down;
	RVec3 point_pos;
	bool move_point = false;
	MatX3 mesh_3d_normals_down;

	Vec2 projected_xh;
	double xh_radius = 0.25;

	bool menus_are_visible = true;
	bool all_menus_are_visible = true;
	Vec2i old_pos_bar, old_viewer_bar_pos, mini_menu_old_pos;

	bool mesh_loaded = false;
	bool mouse_on_uv_side = false;
	bool mouse_over_3d_mesh = false;
	bool mouse_over_2d_mesh = false;
	bool mouse_over_uv_mesh = false;
	bool update_highlights = true;
	bool release_when_translation_done = false;

	// minimalistic menu
	nanogui::Window* mini_menu;
	nanogui::Window* bar_window;

	// Adjacency matrix
	vector<vector<int>> adjacency_list;
	vector<vector<int>> adjacency_list_quad;
	// threshold for separation edge highlighting
	double sep_thresh = 0.1;

	// painting weights
	MatX3 painting_weights, paintint_ones, painting_colors;

	// coords used in process_mouse_move
	int curr_mouse_x, curr_mouse_y;
	bool mouse_updated = false;

	// during coloring just add all newly colored faces into a list, which is later processed
	// to only contain unique elements and then sent to the positional constraint energy
	// this sending and postprocessing is done when releasing the mouse
	// also the position fetching is done then
	vector<int> new_face_hits_fixed, new_face_hits_active, new_face_hits_remove;

	// updated positions of the active triangles
	map<int, Mat32> updated_active_triangles_positions;

	// store the position of the active triangles when coloring them as active for proper offset calculation/animation
	map<int, Mat32> active_triangles_down_positions;

	// maximal distance between each vertex pair in the source mesh
	// needed for ray-mesh intersection to get a decent approximation of
	// the size of the object. this is used to find the value
	// of this variable
	double intersection_delta = 0.;

	// some mouse states
	bool MOUSE_LEFT = false;
	bool MOUSE_MID = false;
	bool MOUSE_RIGHT = false;

	// z position of the crosshair (such that its "above" all vertices)
	double crosshair_z_pos;

	// bounding box stuff
	RVec3 dot_on_mouse_pos;
	vector<pair<RVec2, RVec2>> bbox_corners;
	vector<pair<RVec2, RVec2>> down_corners;
	bool bbox_exists = false;
	bool translate_bbox_corner_1 = false, translate_bbox_corner_2 = false;
	int translate_box_nr = -1;
	bool is_adding_first_corner = true;

	// for easy drawing set these as well (same index for same box)
	vector<RVec3> top_left, bottom_left, bottom_right, top_right;
	double max_perimeter_to_area_ratio=10000;
	double max_distortion=3;
	bool colorByRGB=false;
	MatX3 RGBColors;
	vector<uint> nonOverlappingTriIDs;
	bool largest_overlapping_part_enabled;
	set<pair<int, int>> patchBoundaryEdges;
};

#endif