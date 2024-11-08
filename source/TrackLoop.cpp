// C++ dependencies
#include <iostream>

// X17 dependencies
#include "Points.h"
#include "Reconstruction.h"
#include "TrackLoop.h"
#include "X17Utilities.h"

namespace X17
{
    //// Public methods.

    void TrackLoop::AddTask(RecoTask* task)
    {
        task->m_loop = this;
        m_tasks.push_back(task);
    }

    void TrackLoop::ProcessSingle(TTree* single_track)
    {
        curr_loop = SINGLE;

        curr_micro_tree = single_track;

        curr_micro.SetTTreeBranches(single_track);
        for (RecoTask* t : m_tasks) t->PreElectronLoop();

        // Looping through all ionization electrons.
        int n_electrons = 0;
        for (int i = 0; i < single_track->GetEntries(); i++)
        {
            single_track->GetEntry(i);
            if (IsInSector(curr_micro.x1(),curr_micro.y1(),0) && IsInSector(curr_micro.GetInitPos(),-0.01))
            {
                n_electrons++;
                curr_reco = Reconstruct(*map,curr_micro);
                for (RecoTask* t : m_tasks) t->ElectronLoop();
            }
        }

        std::cout << "\nNumber of electrons in the TPC region: " << n_electrons << "\n";

        for (RecoTask* t : m_tasks) t->PostElectronLoop();
    }

    void TrackLoop::ProcessMulti(TTree* micro_tracks, int n_process)
    {
        curr_loop = MULTI;

        micro_tracks->SetBranchAddress("track_small",&curr_microtrack);
        for (RecoTask* t : m_tasks) t->PreTrackLoop();

        int n_tracks = micro_tracks->GetEntries();
        for (int i = 0; i < n_tracks; i++)
        {
            curr_track_index = i;

            micro_tracks->GetEntry(i);
            // if(curr_microtrack->electron) continue; // ONLY FOR TEST!!!
            if((100 * i) % n_tracks == 0) std::cout << "Progress: " << 100 * i / n_tracks << " \%\n";
            std::cout << "Track " << i+1 << " out of " << n_tracks << ".\n";

            for (RecoTask* t : m_tasks) t->PreElectronLoop();

            int n_electrons = 0;
            for (MicroPoint p : curr_microtrack->points) 
            {
                curr_micro = p;

                // Use only electrons that started and ended in the sector, endpoint not further than 0.5 cm from readout
                if (IsInSector(curr_micro.x1(), curr_micro.y1(), 0) && IsInSector(curr_micro.GetInitPos(), -0.01) && curr_micro.z1() > 7.5)
                {
                    n_electrons++;
                    curr_reco = Reconstruct(*map,curr_micro);
                    for (RecoTask* t : m_tasks) t->ElectronLoop();
                }
            }

            for (RecoTask* t : m_tasks) t->PostElectronLoop();

            if (i == n_process - 1) break;
        }
        
        for (RecoTask* t : m_tasks) t->PostTrackLoop();
    }

    void TrackLoop::ProcessRK(TTree* rk_tracks, int n_process)
    {
        curr_loop = RK;

        rk_tracks->SetBranchAddress("track",&curr_rk);
        for (RecoTask* t : m_tasks) t->PreTrackLoop();

        int n_tracks = rk_tracks->GetEntries();
        for (int i = 0; i < n_tracks; i++)
        {
            rk_tracks->GetEntry(i);
            if((100 * i) % n_tracks == 0) std::cout << "Progress: " << 100 * i / n_tracks << " \%\n";
            std::cout << "Track " << i+1 << " out of " << n_tracks << ".\n";

            for (RecoTask* t : m_tasks) t->PreElectronLoop();

            for (RKPoint p : curr_rk->points) 
            {
                curr_rkpoint = p;
                for (RecoTask* t : m_tasks) t->ElectronLoop();
            }

            for (RecoTask* t : m_tasks) t->PostElectronLoop();

            if (i == n_process - 1) break;
        }
        
        for (RecoTask* t : m_tasks) t->PostTrackLoop();
    }
} // namespace X17
