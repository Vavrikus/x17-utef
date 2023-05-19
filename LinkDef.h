#include "Field.h"
#include "Points.h"

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
    #endif
} // namespace X17