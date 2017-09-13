#include "ApplicationPipeline.h"

#include <igl/viewer/Viewer.h>

using namespace std;
using namespace igl::viewer;

int main(int argc, char** argv)
{
	unsigned int num_threads = max(atoi(getenv("OMP_NUM_THREADS")), 1);
	omp_set_num_threads(num_threads);

	Viewer viewer;
	viewer.callback_init = [&](Viewer& v) {
		// parse command line arguments
		bool is_loaded = false;
		for (int i = 1; i < argc; ++i)
		{
			if (strcmp(argv[i], "-s") == 0)
			{
				cout << "load scene file: " << argv[i + 1] << endl;
				is_loaded = v.load_scene(argv[i + 1]);
				if (is_loaded == false)
				{
					cout << "file not found: " << argv[i + 1] << endl;
				}
			}
		}
		return false;
	};

	// add pipieline plugin
	ApplicationPipeline pipeline;
	viewer.core.rotation_type = igl::viewer::ViewerCore::ROTATION_TYPE_TRACKBALL;
	viewer.plugins.push_back(&pipeline);
	viewer.core.background_color << 1., 1., 1., 1.; // 0.25, 0.25, 0.25, 1.0;
	pipeline.SetArguments(argc, argv);

	// set to ortographic view to synchronize pointer with models
	viewer.core.orthographic = true;

	// start viewer
	viewer.launch(true, false, 2400, 1350);

	return 0;
}