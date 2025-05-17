#pragma once

// C++ dependencies
#include <vector>

// ROOT dependencies
#include "TTree.h"

// X17 dependencies
#include "Track.h"

namespace X17
{
    class TrackLoop;

    /// @brief An abstract class for tasks for TrackLoop.
    class RecoTask
    {
        friend class TrackLoop;
    protected:
        TrackLoop* m_loop = nullptr; // The loop that will run this task.
    public:
        /// @brief To be run before any loops.
        virtual void PreTrackLoop() { }

        /// @brief To be run inside the track loop before the electron loop.
        virtual void PreElectronLoop() { }

        /// @brief To be run during the electron loop on each track.
        virtual void ElectronLoop() { }

        /// @brief To be run inside the track loop.
        virtual void PostElectronLoop() { }

        /// @brief To be run after all loops.
        virtual void PostTrackLoop() { }
    };

    class TrackLoop
    {
    public:
        const Field<MapPoint>& map;            // The ionization electron drift map.
        Field<Vector>* magfield;         // The magnetic field simulated data.

        TTree* curr_micro_tree;          // Current tree with microscopic simulation result.
        MicroPoint curr_micro;           // Current microscopic simulation point.
        RecoPoint curr_reco;             // Current reconstructed point.

        TrackRK* curr_rk;                // Current Runge-Kutta simulated track.
        RKPoint curr_rkpoint;            // Current point on the current Runge-Kutta track.

        TrackMicro* curr_microtrack = 0; // Current microscopically simulated track.
        int curr_track_index;            // Index of the current track in the file(s).

        /// @brief Enumeration of the type of the loop currently running.
        enum LoopType { NONE, SINGLE, MULTI, RK };

        LoopType curr_loop = NONE;    // Type of the current loop.

        bool make_track_plots = true; // Should the tracks be ploted?
        
    private:
        std::vector<RecoTask*> m_tasks; // Vector of all tasks to be run.

    public:
        /// @brief Constructor of TrackLoop.
        /// @param map Pointer to the ionization electron drift map.
        /// @param magfield Pointer to magnetic field data.
        TrackLoop(const Field<MapPoint>& map, Field<Vector>* magfield) : map(map), magfield(magfield) { }

        /// @brief Adds a task to TrackLoop task list.
        /// @param task The task to be added.
        void AddTask(RecoTask* task);

        /// @brief Runs all of the tasks for the single track.
        /// @param single_track TTree with the simulated ionization electrons.
        void ProcessSingle(TTree* single_track);

        /// @brief Runs all of the tasks for tracks simulated by microscopic simulation.
        /// @param micro_tracks TTree with the simulated tracks.
        /// @param n_process The number of tracks to process. Default value is -1, which processes all tracks.
        void ProcessMulti(TTree* micro_tracks, int n_process = -1);

        /// @brief Runs all of the tasks for tracks simulated by Runge-Kutta.
        /// @param rk_tracks TTree with the simulated tracks.
        /// @param n_process The number of tracks to process. Default value is -1, which processes all tracks.
        void ProcessRK(TTree* rk_tracks, int n_process = -1);
    };
} // namespace X17
