// C++ dependenciess
#include <algorithm>
#include <vector>

// ROOT dependencies
#include "TMarker3DBox.h"
// X17 dependencies
#include "RecoPoint.h"
#include "X17Utilities.h"

namespace X17
{
  std::vector<TMarker3DBox*> GetDataMarkers(const std::vector<RecoPoint>& data, double zbin_size)
  {
    std::vector<TMarker3DBox*> markers;
    constexpr float max_size = 0.75F;

    // Find maximal count.
    double max_count = 0;
    for (RecoPoint p : data)
      max_count = std::max<double>(p.count, max_count);

    // Create markers.
    for (RecoPoint p : data)
    {
      using namespace constants;

      float x        = static_cast<float>(p.x());
      float y        = static_cast<float>(p.y());
      float z        = static_cast<float>(p.z());
      float rel_size = max_size * static_cast<float>(p.count / max_count);
      float xlen     = rel_size * static_cast<float>(pad_width / 2.0);
      float ylen     = rel_size * static_cast<float>(pad_height / 2.0);
      float zlen     = rel_size * static_cast<float>(zbin_size / 2.0);

      markers.push_back(new TMarker3DBox(x, y, z, xlen, ylen, zlen, 0, 0));
    }

    return markers;
  }
} // namespace X17