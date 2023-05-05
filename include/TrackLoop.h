#pragma once

#include <iostream>
#include <vector>

#include "TTree.h"

#include "Field.h"
#include "Points.h"
#include "Reconstruction.h"
#include "Track.h"
#include "X17Utilities.h"

namespace X17
{
    class TrackLoop;

    /// @brief An abstract class for tasks for TrackLoop.
    class RecoTask
    {
        friend class TrackLoop;
    protected:
        TrackLoop* loop = nullptr; // The loop that will run this task.
    public:
        /// @brief To be run before any loops.
        virtual void PreTrackLoop()     { }

        /// @brief To be run inside the track loop before the electron loop.
        virtual void PreElectronLoop()  { }

        /// @brief To be run during the electron loop on each track.
        virtual void ElectronLoop()     { }

        /// @brief To be run inside the track loop.
        virtual void PostElectronLoop() { }

        /// @brief To be run after all loops.
        virtual void PostTrackLoop()    { }
    };

    class TrackLoop
    {
    public:
        Field<MapPoint>* map;    // The ionization electron drift map.
        Field<Vector>* magfield; // The magnetic field simulated data.

        TTree* curr_micro_tree;  // Current tree with microscopic simulation result.
        MicroPoint curr_micro;   // Current microscopic simulation point.
        RecoPoint curr_reco;     // Current reconstructed point.

        TrackRK* curr_rk;        // Current Runge-Kutta simulated track.
        RKPoint curr_rkpoint;    // Current point on the current Runge-Kutta track.

        /// @brief Constructor of TrackLoop.
        /// @param map Pointer to the ionization electron drift map.
        /// @param magfield Pointer to magnetic field data.
        TrackLoop(Field<MapPoint>* map, Field<Vector>* magfield) : map(map), magfield(magfield) { }

        /// @brief Adds a task to TrackLoop task list.
        /// @param task The task to be added.
        void AddTask(RecoTask* task)
        {
            task->loop = this;
            tasks.push_back(task);
        }

        /// @brief Runs all of the tasks for the single track.
        /// @param single_track TTree with the simulated ionization electrons.
        void ProcessSingle(TTree* single_track)
        {
            curr_micro_tree = single_track;

            single_track->SetBranchAddress("point",&curr_micro);
            for (RecoTask* t : tasks) t->PreElectronLoop();

            // Looping through all ionization electrons.
            int n_electrons = 0;
            for (int i = 0; i < single_track->GetEntries(); i++)
            {
                single_track->GetEntry(i);
                if (IsInSector(curr_micro.x1,curr_micro.y1,0) && IsInSector(curr_micro.GetInitPos(),-0.01))
                {
                    n_electrons++;
                    curr_reco = Reconstruct(*map,curr_micro);
                    for (RecoTask* t : tasks) t->ElectronLoop();
                }
            }

            std::cout << "\nNumber of electrons in the TPC region: " << n_electrons << "\n";

            for (RecoTask* t : tasks) t->PostElectronLoop();
        }

        /// @brief Runs all of the tasks for tracks simulated by Runge-Kutta.
        /// @param rk_tracks TTree with the simulated tracks.
        /// @param n_process The number of tracks to process. Default value is -1, which processes all tracks.
        void ProcessRK(TTree* rk_tracks, int n_process = -1)
        {
            rk_tracks->SetBranchAddress("track",&curr_rk);
            for (RecoTask* t : tasks) t->PreTrackLoop();

            int n_tracks = rk_tracks->GetEntries();
            for (int i = 0; i < n_tracks; i++)
            {
                rk_tracks->GetEntry(i);
                if((100 * i) % n_tracks == 0) std::cout << "Progress: " << 100*i/n_tracks << " \%\n";
                std::cout << "Track " << i+1 << " out of " << n_tracks << ".\n";

                for (RecoTask* t : tasks) t->PreElectronLoop();

                for (RKPoint p : curr_rk->points) 
                {
                    curr_rkpoint = p;
                    for (RecoTask* t : tasks) t->ElectronLoop();
                }

                for (RecoTask* t : tasks) t->PostElectronLoop();

                if (i == n_process - 1) break;
            }
            
            for (RecoTask* t : tasks) t->PostTrackLoop();
        }
    private:
        std::vector<RecoTask*> tasks; // Vector of all tasks to be run.
    };
} // namespace X17
