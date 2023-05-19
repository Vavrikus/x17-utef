#include "Field.h"
#include "Points.h"
#include "Track.h"

namespace X17
{
    template class X17::Field<X17::MapPoint>;

    #ifdef __CLING__
    #pragma link C++ nestedclasses;
    #pragma link C++ nestedtypedefs;
    #pragma link C++ class X17::Field<X17::MapPoint>+;
    #pragma link C++ class X17::MapPoint+;
    #pragma link C++ class X17::EndPoint+;
    #pragma link C++ class X17::Vector+;
    #pragma link C++ class X17::DataPoint+;
    #pragma link C++ class X17::TrackRK+;
    #pragma link C++ class X17::StartPoint+;
    #endif
} // namespace X17