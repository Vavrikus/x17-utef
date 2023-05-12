// X17 dependencies
#include "TrackLoop.h"

namespace X17
{
    void TrackLoop::AddTask(RecoTask* task)
    {
        task->loop = this;
        tasks.push_back(task);
    }

    void TrackLoop::ProcessSingle(TTree* single_track)
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

    void TrackLoop::ProcessRK(TTree* rk_tracks, int n_process)
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
} // namespace X17
