#pragma once

#ifndef APPLICATION_PIPELINE_STEP_H
#define APPLICATION_PIPELINE_STEP_H

#include <igl/viewer/Viewer.h>
#include <nanogui/formhelper.h>
#include <nanogui/window.h>

#include "ApplicationPipeline.h"

struct DataObject;
namespace igl { namespace viewer { class Viewer; } }
namespace nanogui { class FormScreen; }

class ApplicationPipelineStep
{
  friend class ApplicationPipeline;
public:

  // intialization
  void Init(igl::viewer::Viewer *viewer) { this->viewer = viewer; }

  // viewer shuts down
  virtual void Shutdown() {}

  // saving and loading of a mesh
  virtual bool Load(std::string filename) { return false; }
  virtual bool Save(std::string filename) { return false; }
  virtual bool PostLoad() { return false; }

  // de/serialization events
  virtual bool Serialize(std::vector<char>& buffer) const { return false; }
  virtual bool Deserialize(const std::vector<char>& buffer) { return false; }

  // draw stuff
  virtual bool PreDraw() { return false; }
  virtual bool PostDraw() { return false; }

  // mouse events
  virtual bool MouseDown(int button,int modifier) { return false; }
  virtual bool MouseUp(int button,int modifier) { return false; }
  virtual bool MouseMove(int mouse_x,int mouse_y) { return false; }
  virtual bool MouseScroll(float delta_y) { return false; }

  // keyboard events
  virtual bool KeyDown(int key,int modifiers) { return false; }
  virtual bool KeyUp(int key,int modifiers) { return false; }
  virtual bool KeyPressed(int key,int modifiers) { return false; }

  // set and init first time initialization stuff (add menu entries)
  virtual void InitPipelineStep() {}

  // set and init pipeline state depending stuff (state switch)
  virtual void EnterPipelineStep() {}

  // set and init pipeline state depending stuff (state switch)
  virtual void LeavePipelineStep() {}

  // compute the next step of the algorithm!
  virtual bool ComputePipelineStep() { return false; }

protected:

  // access to the igl viewer
  igl::viewer::Viewer *viewer;

  // access to application pipeline
  ApplicationPipeline* pipeline;

  // access to the menu
  nanogui::FormHelper* bar;
	nanogui::Window* orig_window;

  // access to the data object
  DataObject* data;

};

// defines serialization for ApplicationPipelineStep
namespace igl
{
  namespace serialization
  {
    template<> inline void serialize(const ApplicationPipelineStep& obj,std::vector<char>& buffer)
    {
      obj.Serialize(buffer);
    }

    template<> inline void deserialize(ApplicationPipelineStep& obj,const std::vector<char>& buffer)
    {
      obj.Deserialize(buffer);
    }
  }
}

#endif
