// X17 dependencies
#include "RecoPoint.h"
#include "X17Utilities.h"


namespace X17
{
    std::vector<TMarker3DBox*> GetDataMarkers(const std::vector<RecoPoint> &data, double zbin_size)
    {
        std::vector<TMarker3DBox*> markers;
        constexpr double max_size = 0.75;

        // Find maximal count.
        int max_count = 0;
        for (RecoPoint p : data) if (p.count > max_count) max_count = p.count;

        // Create markers.
        for (RecoPoint p : data)
        {
            using namespace constants;

            double rel_size = max_size * p.count / max_count;
            double xlen = rel_size * pad_width  / 2.0;
            double ylen = rel_size * pad_height / 2.0;
            double zlen = rel_size * zbin_size  / 2.0;

            markers.push_back(new TMarker3DBox(p.x(),p.y(),p.z(),xlen,ylen,zlen,0,0));
        }
        
        return markers;
    }
} // namespace X17