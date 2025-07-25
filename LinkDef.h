#include "Field.h"
#include "Points.h"
#include "Track.h"

namespace X17
{
    template class X17::Field<X17::MapPoint>;
    template class X17::Matrix<4,4>;

    #pragma link C++ nestedclasses;
    #pragma link C++ nestedtypedefs;
    #pragma link C++ class X17::Vector+;
    #pragma link C++ class X17::Field<X17::Vector>+;
    #pragma link C++ class X17::StartPoint+;
    #pragma link C++ class X17::EndPoint+;
    #pragma link C++ class X17::EndPointDiscrete+;
    #pragma link C++ class X17::MicroPoint+;
    #pragma link C++ class X17::RKPoint+;
    #pragma link C++ class std::vector<X17::RKPoint>+;
    #pragma link C++ class X17::MapPoint+;
    #pragma link C++ class X17::DataPoint+;
    #pragma link C++ class X17::DriftLinePoint+;
    #pragma link C++ class X17::Field<X17::MapPoint>+;
    #pragma link C++ class X17::TrackRK+;
    #pragma link C++ class X17::TrackMicro+;
    #pragma link C++ class X17::TrackInfo+;
    #pragma link C++ class X17::Matrix<4,4>+;
} // namespace X17