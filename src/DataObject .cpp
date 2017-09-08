#include "DataObject.h"
#include <igl/viewer/ViewerData.h>
#include <igl/serialize.h>

DataObject::DataObject()
{
}

SERIALIZE_TYPE_SOURCE(DataObject,
	SERIALIZE_MEMBER_NAME(V, "V")
	SERIALIZE_MEMBER_NAME(F, "F")
)