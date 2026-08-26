#pragma once
// IWYU pragma: private, include "Points.h"

// ROOT dependencies
#include "Rtypes.h"
#include "TChain.h"
#include "TTree.h"

// X17 dependencies
#include "EndPoint.h"
#include "StartPoint.h"

namespace X17
{
  /// @brief A struct for storing the results of microscopic simulation of ionization electrons.
  struct MicroPoint
  {
    StartPoint start; // Initial coordinates of the electron ([cm] and [ns]).
    EndPoint end;     // Final coordinates of the electron ([cm] and [ns]).
    double e0 = 0.0;  // Initial energy [eV].
    double e1 = 0.0;  // Final energy [eV].

    /// @brief The default constructor. Initializes time to -1, everything else to 0.
    MicroPoint()
      : start(), end()
    {
    }

    /// @brief Getter for the initial x variable.
    /// @return Initial x-coordinate [cm].
    [[nodiscard]] double x0() const { return start.x(); }

    /// @brief Getter for the initial y variable.
    /// @return Initial y-coordinate [cm].
    [[nodiscard]] double y0() const { return start.y(); }

    /// @brief Getter for the initial z variable.
    /// @return Initial z-coordinate [cm].
    [[nodiscard]] double z0() const { return start.z(); }

    /// @brief Getter for the initial t variable.
    /// @return Initial time [ns].
    [[nodiscard]] double t0() const { return start.t; }

    /// @brief Getter for the final x variable.
    /// @return Final x-coordinate [cm].
    [[nodiscard]] double x1() const { return end.x(); }

    /// @brief Getter for the final y variable.
    /// @return Final y-coordinate [cm].
    [[nodiscard]] double y1() const { return end.y(); }

    /// @brief Getter for the final z variable.
    /// @return Final z-coordinate [cm].
    [[nodiscard]] double z1() const { return end.z(); }

    /// @brief Getter for the final t variable.
    /// @return Final time [ns].
    [[nodiscard]] double t1() const { return end.t; }

    /// @brief Returns the initial position of the ionization electron.
    /// @return The initial position of the electron.
    [[nodiscard]] Vector GetInitPos() const { return this->start.point; }

    /// @brief Creates the branches of the given tree used for output.
    /// @param tree TTree used for output.
    void MakeTTreeBranches(TTree* tree);

    /// @brief Sets the branches of a TChain to this MicroPoint object.
    /// @param chain Pointer to TChain object.
    /// @param old_data Boolean indicating whether the data is from old simulations with different coordinate system
    /// (zxy). Defaults to true.
    void SetTChainBranches(TChain* chain, bool old_data = true);

    /// @brief Sets the branches of a TTree to this MicroPoint object. Uses different branch structure than
    /// MakeTTreeBranches.
    /// @param tree Pointer to a TTree object.
    void SetTTreeBranches(TTree* tree);

    ClassDefNV(MicroPoint, 1);
  };
} // namespace X17