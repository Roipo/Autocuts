#pragma once

#ifndef APPLICATION_PIPELINE_GROUP_H
#define APPLICATION_PIPELINE_GROUP_H

#include <igl\viewer\Viewer.h>

struct DataObject;
namespace igl { namespace viewer { class Viewer; } }

class ApplicationPipelineGroup
{
  friend class ApplicationPipeline;
public:

  ApplicationPipelineGroup() : initialized(false) { };

  // viewer shuts down
  virtual void Shutdown() {};

  // saving and loading of a mesh
  virtual bool Load(std::string filename) { return false; };
  virtual bool Save(std::string filename) { return false; };
  virtual bool PostLoad() { return false; };

  // draw stuff
  virtual bool PreDraw() { return false; };
  virtual bool PostDraw() { return false; };

  // mouse events
  virtual bool MouseDown(int button,int modifier) { return false; };
  virtual bool MouseUp(int button,int modifier) { return false; };
  virtual bool MouseMove(int mouse_x,int mouse_y) { return false; };
  virtual bool MouseScroll(float delta_y) { return false; };

  // keyboard events
  virtual bool KeyDown(int key,int modifiers) { return false; };
  virtual bool KeyUp(int key,int modifiers) { return false; };
  virtual bool KeyPressed(int key,int modifiers) { return false; };

  // set and init first time initialization stuff (add menu entries)
  virtual void InitPipelineGroup() {};

  // set and init pipeline state depending stuff (state switch)
  virtual void PipelineGroupChanged() {};

protected:

  // access to the igl viewer
  igl::viewer::Viewer *viewer;

  // access to application pipeline
  ApplicationPipeline* pipeline;

  // access to the menu
  nanogui::FormHelper* bar;

  // acces to the data object
  DataObject* data;

private:

  bool initialized;

};


#endif