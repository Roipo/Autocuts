#include "ApplicationPipeline.h"
#include "ApplicationPipelineStep.h"
#include "ApplicationPipelineGroup.h"
#include "DataObject.h"

#include <igl/viewer/Viewer.h>
#include <igl/viewer/ViewerData.h>
#include <igl/serialize.h>
#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>

using namespace std;

// initialization (runs every time a mesh is loaded or cleared)
void ApplicationPipeline::init(igl::viewer::Viewer *viewer)
{
  if(bar == nullptr)
  {
    ViewerPlugin::init(viewer);

    initSteps();
    initMenu();

    // init pipeline step groups
    for(auto&& v : pipelineGroups)
    {
      if(v != nullptr && v->initialized == false)
      {
        v->viewer = viewer;
        v->pipeline = this;
        v->bar = bar;
        v->InitPipelineGroup();
        v->initialized = true;
      }
    }

    // init pipeline steps
    pipelineData.resize(pipelineSteps.size());
    pipelineViewerData.resize(pipelineSteps.size());
    for(int i=0;i<pipelineSteps.size();i++)
    {
      pipelineSteps[i]->viewer = viewer;
      pipelineSteps[i]->pipeline = this;
      pipelineSteps[i]->bar = bar;
			pipelineSteps[i]->orig_window = orig_window;

      pipelineData[i] = new DataObject();
      pipelineViewerData[i] = new std::vector<igl::viewer::ViewerData>();
      pipelineSteps[i]->data = pipelineData[i];
      
      pipelineSteps[i]->InitPipelineStep();
    }

    // init pipeline buffers
    pipelineDataBufferStepId.clear();
    pipelineViewerDataBuffer.clear();
    pipelineDataBuffer.clear();
    for(int i=0;i<4;i++)
    {
      pipelineDataBufferStepId.push_back(-1);
      pipelineDataBuffer.push_back(new DataObject());
      pipelineViewerDataBuffer.push_back(new std::vector<igl::viewer::ViewerData>());
    }

    if(pipelineData.size() > 0)
    {
      data = pipelineData[stepId];
      if(pipelineGroups[stepId] != nullptr)
        pipelineGroups[stepId]->data = data;
      pipelineSteps[stepId]->EnterPipelineStep();
    }

    // Parse command line arguments
    glfwPollEvents();
    if(glfwGetKey(viewer->window,GLFW_KEY_LEFT_SHIFT) != GLFW_PRESS)
    {
      bool isLoaded = false;
      for(int i=1;i<argc;i++)
      {
        if(strcmp(argv[i],"-p")==0) {
          cout << "load pipeline state file: " << argv[i+1] << endl;
          isLoaded = LoadPipelineState(argv[i+1]);
          if(isLoaded == false)
            cout << "file not found: " << argv[i+1] << endl;
        }
      }
    }
  }
}

// called before shutdown of viewer
void ApplicationPipeline::shutdown()
{
}

// called when the scene is serialized
bool ApplicationPipeline::serialize(std::vector<char>& buffer) const
{
  bool serialized = false;

  for(auto&& v : pipelineSteps)
  {
    serialized |= v->Serialize(buffer);
  }

  for(auto&& v : pipelineGroups)
  {
  }

  return serialized;
}

// called when the scene is deserialized
bool ApplicationPipeline::deserialize(const std::vector<char>& buffer)
{
  bool deserialized = false;

  for(auto&& v : pipelineSteps)
  {
    deserialized |= v->Deserialize(buffer);
  }

  for(auto&& v : pipelineGroups)
  {
  }

  return deserialized;
}

// called before a mesh is loaded
bool ApplicationPipeline::load(std::string filename)
{
  if(pipelineGroups.size() > 0 && pipelineGroups[stepId] != nullptr)
    pipelineGroups[stepId]->Load(filename);

  if(pipelineSteps.size() > 0)
    return pipelineSteps[stepId]->Load(filename);

  return false;
}

// called before a mesh is saved
bool ApplicationPipeline::save(std::string filename)
{
  if(pipelineGroups.size() > 0 && pipelineGroups[stepId] != nullptr)
    pipelineGroups[stepId]->Save(filename);

  if(pipelineSteps.size() > 0)
    return pipelineSteps[stepId]->Save(filename);

  return false;
}

// runs immediately after a new mesh has been loaded.
bool ApplicationPipeline::post_load()
{
  if(pipelineGroups.size() > 0 && pipelineGroups[stepId] != nullptr)
    pipelineGroups[stepId]->PostLoad();

  if(pipelineSteps.size() > 0)
    return pipelineSteps[stepId]->PostLoad();

  return false;
}

// called before the draw procedure of viewer
bool ApplicationPipeline::pre_draw()
{
  stepIdGUI = stepId+1;

  if(pipelineGroups.size() > 0 && pipelineGroups[stepId] != nullptr)
    pipelineGroups[stepId]->PreDraw();
  
  if(pipelineSteps.size() > 0)
    return pipelineSteps[stepId]->PreDraw();

  return false;
}

// called after the draw procedure of viewer
bool ApplicationPipeline::post_draw()
{
  if(pipelineGroups.size() > 0 && pipelineGroups[stepId] != nullptr)
    pipelineGroups[stepId]->PostDraw();
  
  if(pipelineSteps.size() > 0)
    return pipelineSteps[stepId]->PostDraw();

  return false;
}

// called when the mouse button is pressed
bool ApplicationPipeline::mouse_down(int button,int modifier)
{
  bool handled = false;
  
  if(pipelineGroups.size() > 0 && pipelineGroups[stepId] != nullptr)
    handled = pipelineGroups[stepId]->MouseDown(button,modifier);
  
  if(pipelineSteps.size() > 0)
    handled = handled || pipelineSteps[stepId]->MouseDown(button,modifier);
  
  return handled;
}

// called when the mouse button is released
bool ApplicationPipeline::mouse_up(int button,int modifier)
{
  bool handled = false;
  
  if(pipelineGroups.size() > 0 && pipelineGroups[stepId] != nullptr)
    handled = pipelineGroups[stepId]->MouseUp(button,modifier);
  
  if(pipelineSteps.size() > 0)
    handled = handled || pipelineSteps[stepId]->MouseUp(button,modifier);
  
  return handled;
}

// called every time the mouse is moved
// - mouse_x and mouse_y are the new coordinates of the mouse pointer in screen coordinates
bool ApplicationPipeline::mouse_move(int mouse_x,int mouse_y)
{
  bool handled = false;
  
  if(pipelineGroups.size() > 0 && pipelineGroups[stepId] != nullptr)
    handled = pipelineGroups[stepId]->MouseMove(mouse_x,mouse_y);
  
  if(pipelineSteps.size() > 0)
    handled = handled || pipelineSteps[stepId]->MouseMove(mouse_x,mouse_y);
  
  return handled;
}

// called every time the scroll wheel is moved
bool ApplicationPipeline::mouse_scroll(float delta_y)
{
  bool handled = false;
  
  if(pipelineGroups.size() > 0 && pipelineGroups[stepId] != nullptr)
    handled = pipelineGroups[stepId]->MouseScroll(delta_y);
  
  if(pipelineSteps.size() > 0)
    handled = handled || pipelineSteps[stepId]->MouseScroll(delta_y);
  
  return handled;
}

// called when a keyboard key is pressed
bool ApplicationPipeline::key_pressed(unsigned int key,int modifiers)
{
  bool handled = false;

  if(pipelineGroups.size() > 0 && pipelineGroups[stepId] != nullptr)
    handled = pipelineGroups[stepId]->KeyPressed(key,modifiers);

  if(pipelineSteps.size() > 0)
    handled = handled || pipelineSteps[stepId]->KeyPressed(key,modifiers);

  return false;
}

// called when a keyboard key is pressed
bool ApplicationPipeline::key_down(int key,int modifiers)
{
  int keyStateId = key - GLFW_KEY_F1;

  // Check for Fxx keys
  if(key >= GLFW_KEY_F1 && key <= GLFW_KEY_F8)
  {
    if(modifiers & GLFW_MOD_SHIFT)
    {
      SwitchPipelineState(keyStateId);
    }
    else
    {
      computeStep(keyStateId);
    }
    return true;
  }
  else if(key >= GLFW_KEY_F9 && key <= GLFW_KEY_F12)
  {
    if(modifiers & GLFW_MOD_SHIFT)
    {
      // copy pipeline state to buffer
      BufferPipelineState(keyStateId-8);
    }
    else
    {
      // load pipeline state from buffer
      LoadBufferedPipelineState(keyStateId-8);
    }
    return true;
  }

  bool handled = false;
  
  if(pipelineGroups.size() > 0 && pipelineGroups[stepId] != nullptr)
    handled = pipelineGroups[stepId]->KeyDown(key,modifiers);
  
  if(pipelineSteps.size() > 0)
    handled = handled || pipelineSteps[stepId]->KeyDown(key,modifiers);
  
  return handled;
}

// called when a keyboard key is release
bool ApplicationPipeline::key_up(int key,int modifiers)
{
  bool handled = false;

  if(pipelineGroups.size() > 0 && pipelineGroups[stepId] != nullptr)
    handled = pipelineGroups[stepId]->KeyUp(key,modifiers);

  if(pipelineSteps.size() > 0)
    handled = handled || pipelineSteps[stepId]->KeyUp(key,modifiers);

  return handled;
}
  
ApplicationPipeline::ApplicationPipeline()
{ 
  plugin_name = "Pipeline";

  stepId = 0;
  stepIdGUI = stepId;

  bar = nullptr;
  data = nullptr;

  isSaveViewerSettings = true;
}

ApplicationPipeline::~ApplicationPipeline()
{
}

void ApplicationPipeline::SetArguments(int argc,char** argv)
{
  this->argc = argc;
  this->argv = argv;
}

void ApplicationPipeline::initMenu()
{
  using namespace nanogui;

  // init menu bar
  if(bar == nullptr)
  {
    // Create a new nanogui window
    bar = viewer->ngui;
    
    //bar->setInputCellSize(Eigen::Vector2i(80,20));
		orig_window = bar->window();
    Window* window = bar->addWindow(Eigen::Vector2i(10,10),"Application Pipeline");

    bar->addVariable("State",stepIdGUI,false);

    bar->addGroup("Pipeline State");

    Widget* container = new Widget(window);
    bar->addWidget("",container);

    GridLayout* layout = new GridLayout(Orientation::Horizontal,2,Alignment::Fill);
    container->setLayout(layout);

    Button* loadButton = new Button(container,"Load");
    loadButton->setFixedHeight(25);
    loadButton->setCallback([&](){this->LoadPipelineState();});

    Button* saveButton = new Button(container,"Save");
    saveButton->setFixedHeight(25);
    saveButton->setCallback([&](){this->SaveFullPipelineState();});

    bar->addVariable("Include Menu Settings",isSaveViewerSettings,true);

    viewer->screen->performLayout();
  }
}

bool ApplicationPipeline::computeStep(int newStepId)
{
  if(newStepId >= pipelineSteps.size())
  {
    cout << "Pipeline: step " << newStepId << " is not available!" << endl;
    return false;
  }
  
  if(stepId <= newStepId)
  {
    // compute all steps up to newStepId
    for(int s=stepId;s<newStepId;s++)
    {
      if(pipelineData[s] != nullptr)
      {
        // deinit prev step
        pipelineSteps[stepId]->LeavePipelineStep();

        stepId = s+1;

        // load and copy prev state
        if(pipelineData[stepId] == nullptr)
        {
          pipelineData[stepId] = new DataObject();
          pipelineSteps[stepId]->data = pipelineData[stepId];
          if(pipelineGroups[stepId] != nullptr)
            pipelineGroups[stepId]->data = pipelineData[stepId];
        }
        *pipelineData[stepId] = *pipelineData[s];
        data = pipelineData[stepId];
        if(pipelineGroups[stepId] != nullptr)
          pipelineGroups[stepId]->data = data;

        // init next step
        pipelineSteps[stepId]->EnterPipelineStep();

        // compute next step
        pipelineSteps[stepId]->ComputePipelineStep();
      }
      else
        cout << "Pipeline: state " << s << " has not been computed!" << endl;
    }
  }
  else if(newStepId > 0)
  {
    if(pipelineData[newStepId] != nullptr)
    {
      // deinit prev step
      pipelineSteps[stepId]->LeavePipelineStep();
      
      stepId = newStepId;

      // load and copy prev state
      *pipelineData[stepId] = *pipelineData[stepId-1];
      data = pipelineData[stepId];
      if(pipelineGroups[stepId] != nullptr)
        pipelineGroups[stepId]->data = data;

      // init next step
      pipelineSteps[stepId]->EnterPipelineStep();

      // recompute provided step
      pipelineSteps[stepId]->ComputePipelineStep();
    }
    else
      cout << "Pipeline: state " << newStepId << " has not been computed!" << endl;
  }
  else
  {
    SwitchPipelineState(newStepId);
  }

  return true;
}

int ApplicationPipeline::GetPipelineState()
{
  return stepId;
}

void ApplicationPipeline::SwitchPipelineState(int newStateId)
{
  if(pipelineData[newStateId] != nullptr)
  {
    pipelineSteps[stepId]->LeavePipelineStep();
    stepId = newStateId;
    data = pipelineData[stepId];
    if(pipelineGroups[stepId] != nullptr)
      pipelineGroups[stepId]->data = data;
    pipelineSteps[stepId]->data = data;
    pipelineSteps[stepId]->EnterPipelineStep();
  }
  else
  {
    cout << "Pipeline: state " << newStateId << " is not yet computed!" << endl;
  }
}

void ApplicationPipeline::BufferPipelineState(int bufferId)
{
  viewer->data_buffer[viewer->active_data_id] = viewer->data;

  pipelineDataBufferStepId[bufferId] = stepId;
  *pipelineDataBuffer[bufferId] = *data;
  *pipelineViewerDataBuffer[bufferId] = viewer->data_buffer;
}

void ApplicationPipeline::LoadBufferedPipelineState(int bufferId)
{
  if(pipelineDataBufferStepId[bufferId] > -1)
  {
    pipelineSteps[stepId]->LeavePipelineStep();
    
    stepId = pipelineDataBufferStepId[bufferId];
    *data = *pipelineDataBuffer[bufferId];
    if(pipelineGroups[stepId] != nullptr)
      pipelineGroups[stepId]->data = data;
    *pipelineSteps[stepId]->data = *data;
    
    viewer->data_buffer = *pipelineViewerDataBuffer[bufferId];
    viewer->data = viewer->data_buffer[viewer->active_data_id];
    for(auto& v: viewer->data_buffer)
      v.dirty = v.DIRTY_ALL;
    viewer->data.dirty = viewer->data.DIRTY_ALL;

    pipelineSteps[stepId]->EnterPipelineStep();
  }
  else
  {
    cout << "Pipeline: buffer " << bufferId << " is empty!" << endl;
  }
}

void ApplicationPipeline::PrevPipelineState()
{
  if(stepId >= 1)
  {
    pipelineSteps[stepId]->LeavePipelineStep();
    stepId--;
    data = pipelineData[stepId];
    if(pipelineGroups[stepId] != nullptr)
      pipelineGroups[stepId]->data = data;
    pipelineSteps[stepId]->EnterPipelineStep();
  }
}

void ApplicationPipeline::NextPipelineState()
{
  if(stepId < pipelineData.size()-1)
  {
    pipelineSteps[stepId]->LeavePipelineStep();
    stepId++;
    data = pipelineData[stepId];
    if(pipelineGroups[stepId] != nullptr)
      pipelineGroups[stepId]->data = data;
    pipelineSteps[stepId]->EnterPipelineStep();
  }
}

bool ApplicationPipeline::LoadPipelineState()
{
  std::string fname = igl::file_dialog_open();
  return LoadPipelineState(fname.c_str());
}

bool ApplicationPipeline::LoadPipelineState(const std::string& fileName)
{
  const char* fname = fileName.c_str();

  if(fname[0] == 0)
    return false;

  igl::deserialize(stepId,"PipelineStateId",fname);
  igl::deserialize(pipelineData,"PipelineData",fname);
  igl::deserialize(viewer->active_data_id,"ActiveDataId",fname);
  igl::deserialize(viewer->data_buffer,"ViewerData",fname);
  igl::deserialize(viewer->data_ids,"ViewerDataIds",fname);

  viewer->data = viewer->data_buffer[viewer->active_data_id];
  
  // add missing opengl states
  for(int i=viewer->opengl.size();i<viewer->data_buffer.size();i++)
  {
    viewer->opengl.push_back(igl::viewer::OpenGL_state());
    viewer->opengl[viewer->opengl.size()-1].init();
  }

  // remove remaining opengl states
  for(int i=viewer->data_buffer.size();i<viewer->opengl.size();i++)
  {
    viewer->opengl[i].free();
    viewer->opengl.erase(viewer->opengl.begin()+i);
  }

  // init mesh combobox
  if(viewer->currentDataCB != nullptr)
  {
    viewer->currentDataCB->setItems(viewer->data_ids);
  }

  if(isSaveViewerSettings)
    igl::deserialize(pipelineSteps,"PipelineSteps",fname);
  
  if(pipelineData.size() < pipelineSteps.size())
  {
    // add missing states
    vector<DataObject*> tempData = pipelineData;
    pipelineData.resize(pipelineSteps.size(),nullptr);
    for(int i=0;i<tempData.size();i++)
      pipelineData[i] = tempData[i];
  }
    
  for(int i=0;i<pipelineSteps.size();i++)
  {
      pipelineSteps[i]->data = pipelineData[i];
      if(pipelineGroups[i] != nullptr)
        pipelineGroups[i]->data = pipelineData[i];
  }
    
  SwitchPipelineState(stepId);
  
  return true;
}

bool ApplicationPipeline::SaveCurrentPipelineState()
{
  std::string fname = igl::file_dialog_save();
  return SaveCurrentPipelineState(fname);
}

bool ApplicationPipeline::SaveCurrentPipelineState(const std::string& fname)
{
  std::vector<DataObject*> saveData = pipelineData;

  for(int i=0;i<saveData.size();i++)
  {
    if(i != stepId)
      saveData[i] = nullptr;
  }

  return savePipelineState(saveData,fname);
}

bool ApplicationPipeline::SaveFullPipelineState()
{
  std::string fname = igl::file_dialog_save();
  return savePipelineState(pipelineData,fname);
}

bool ApplicationPipeline::SaveFullPipelineState(const std::string& fname)
{
  return savePipelineState(pipelineData,fname);
}

bool ApplicationPipeline::savePipelineState(std::vector<DataObject*>& saveData, const std::string& fileName)
{
  const char* fname = fileName.c_str();

  if(fname[0] == 0)
    return false;

  viewer->data_buffer[viewer->active_data_id] = viewer->data;

  igl::serialize(stepId,"PipelineStateId",fname,true);
  igl::serialize(saveData,"PipelineData",fname);
  igl::serialize(viewer->active_data_id,"ActiveDataId",fname);
  igl::serialize(viewer->data_buffer,"ViewerData",fname);
  igl::serialize(viewer->data_ids,"ViewerDataIds",fname);

  if(isSaveViewerSettings)
    igl::serialize(pipelineSteps,"PipelineSteps",fname);

  return true;
}

void ApplicationPipeline::addStep(ApplicationPipelineStep* step)
{
  addStep(step,nullptr);
}

void ApplicationPipeline::addStep(ApplicationPipelineStep* step, ApplicationPipelineGroup* group)
{
  pipelineSteps.push_back(step);
  pipelineGroups.push_back(group);
}

#include "SolverPlugin.h"

void ApplicationPipeline::initSteps()
{
	addStep(new SolverPlugin());
}