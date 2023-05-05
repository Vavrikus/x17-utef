#pragma once

#include "Field.h"
#include "Points.h"

/// @brief An abstract class representing a map task.
class MapTask
{
protected:
    X17::Field<X17::MapPoint>* map; // A pointer to the ionization electron map.

    /// @brief The constructor of MapTask.
    /// @param map A pointer to the ionization electron drift map.
    MapTask(X17::Field<X17::MapPoint>* map) : map(map) { }

public:
    /// @brief A virtual function to be called before all loops of the map points.
    virtual void PreLoop() { }

    /// @brief A virtual function to be called at the start of the z-loop of the map points.
    /// @param z The z coordinate of the current map point.
    virtual void Z_Loop_Start(const double& z) { }

    /// @brief A virtual function to be called during the xyz-loop of the map points.
    /// @param x The x coordinate of the current map point.
    /// @param y The y coordinate of the current map point.
    /// @param z The z coordinate of the current map point.
    /// @param current The current map point.
    virtual void XYZ_Loop(const double& x, const double& y, const double& z, const X17::MapPoint& current) { }

    /// @brief A virtual function to be called at the end of the z-loop of the map points.
    virtual void Z_Loop_End() { }

    /// @brief A virtual function to be called during the zx-loop of the map points.
    /// @param z The z coordinate of the current map point.
    /// @param x The x coordinate of the current map point.
    /// @param current The current map point.
    virtual void ZX_Loop(const double& z, const double& x, const X17::MapPoint& current) { }

    /// @brief A virtual function to be called after all loops of the map points.
    virtual void PostLoop() { }
};