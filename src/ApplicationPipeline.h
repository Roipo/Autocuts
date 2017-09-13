#pragma once

#ifndef APPLICATION_PIPELINE_H
#define APPLICATION_PIPELINE_H

#include <igl/viewer/Viewer.h>
#include <igl/viewer/ViewerPlugin.h>
#include <nanogui/window.h>
#include <nanogui/formhelper.h>

struct DataObject;
class ApplicationPipelineStep;
class ApplicationPipelineGroup;
class ApplicationPipelineState;
namespace igl { namespace viewer { class ViewerData; } }
namespace nanogui { class FormScreen; }

class ApplicationPipeline : public igl::viewer::ViewerPlugin
{
public:

  // --------------------------------------------------------------------------------
  // implementation of interface of viewer
  // --------------------------------------------------------------------------------

  // initialization
  void init(igl::viewer::Viewer *viewer);

  // called before shutdown of viewer
  void shutdown() override;

  // called before a mesh is loaded
  bool load(std::string filename) override;

  // called before a mesh is saved
  bool save(std::string filename) override;

  // called when the scene is serialized
  bool serialize(std::vector<char>& buffer) const override;

  // called when the scene is deserialized
  bool deserialize(const std::vector<char>& buffer) override;

  // runs immediately after a new mesh has been loaded.
  bool post_load() override;

  // called before the draw procedure of viewer
  bool pre_draw() override;

  // called after the draw procedure of viewer
  bool post_draw() override;

  // called when the mouse button is pressed
  // - button can be GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON or GLUT_RIGHT_BUTTON
  // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
  bool mouse_down(int button,int modifier) override;

  // called when the mouse button is released
  // - button can be GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON or GLUT_RIGHT_BUTTON
  // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
  bool mouse_up(int button,int modifier) override;

  // called every time the mouse is moved
  // - mouse_x and mouse_y are the new coordinates of the mouse pointer in screen coordinates
  bool mouse_move(int mouse_x,int mouse_y) override;

  // called every time the scroll wheel is moved
  // Note: this callback is not working with every glut implementation
  bool mouse_scroll(float delta_y) override;

  // This function is called when a keyboard key is pressed. Unlike key_down
  // this will reveal the actual character being sent (not just the physical
  // key)
  // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
  bool key_pressed(unsigned int key,int modifiers) override;

  // called when a keyboard key is pressed
  // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
  bool key_down(int key,int modifiers) override;

  // called when a keyboard key is release
  // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
  bool key_up(int key,int modifiers) override;

  // --------------------------------------------------------------------------------
  // Application Pipeline
  // --------------------------------------------------------------------------------

  ApplicationPipeline();
  ~ApplicationPipeline();

  void SetArguments(int argc,char** argv);

  // load pipeline state
  bool LoadPipelineState();
  bool LoadPipelineState(const std::string& fname);

  // save pipeline state
  bool SaveCurrentPipelineState();
  bool SaveCurrentPipelineState(const std::string& fname);
  bool SaveFullPipelineState();
  bool SaveFullPipelineState(const std::string& fname);

  // switch between pipeline steps
  int GetPipelineState();
  void SwitchPipelineState(int newStateId);
  void BufferPipelineState(int bufferId);
  void LoadBufferedPipelineState(int bufferId);
  void PrevPipelineState();
  void NextPipelineState();

private:

  // Pointer to the nano gui
	nanogui::FormHelper* bar;
	nanogui::Window* orig_window;

  // Data object contains the whole state of the algorithm
  DataObject* data;

  // arguments
  int argc;
  char** argv;

  std::vector<DataObject*> pipelineData;
  std::vector<std::vector<igl::viewer::ViewerData>*> pipelineViewerData;
  std::vector<DataObject*> pipelineDataBuffer;
  std::vector<std::vector<igl::viewer::ViewerData>*> pipelineViewerDataBuffer;

  int stepId;
  int stepIdGUI; // stepId+1 to fit 'Fxx' keys
  std::vector<int> pipelineDataBufferStepId;
  std::vector<ApplicationPipelineStep*> pipelineSteps;
  std::vector<ApplicationPipelineGroup*> pipelineGroups;

  // GUI flags
  bool isSaveViewerSettings;

  void initMenu();

  bool computeStep(int stepId);
  
  // save pipeline state
  bool savePipelineState(std::vector<DataObject*>& saveData, const std::string& fileName);

  // add a pipeline step
  void addStep(ApplicationPipelineStep* step);
  void addStep(ApplicationPipelineStep* step, ApplicationPipelineGroup* group);

  // here you can add all your algorithm steps
  void initSteps();
};

#endif
