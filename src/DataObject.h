#pragma once

#ifndef DATA_OBJECT_H
#define DATA_OBJECT_H


#include <igl/viewer/Viewer.h>
#include <vector>
#include "igl/serialize.h"
#include "EigenTypes.h"

struct DataObject;
struct Handle;

namespace igl {
  namespace viewer { class ViewerData; }
}

struct DataObject
{
  DataObject();

	MatX3 V;
	MatX3i F;
};

namespace igl {
  namespace serialization
  {
    void _serialize(const DataObject& obj,std::vector<char>& buffer);
    void _deserialize(DataObject& obj,const std::vector<char>& buffer);

    template<> inline void serialize(const DataObject& obj,std::vector<char>& buffer)
    {
      igl::serialization::_serialize(obj,buffer);
    }

    template<> inline void deserialize(DataObject& obj,const std::vector<char>& buffer)
    {
      igl::serialization::_deserialize(obj,buffer);
    }
  }
}

#endif
