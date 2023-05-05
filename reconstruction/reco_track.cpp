#include "TrackLoop.h"

int reco_track()
{
    TrackLoop* loop = new TrackLoop();

    loop->LoadSingle("../data/single_track/new_9010/electrons.root");
    loop->LoadMap("../data/ion_map/map.root");
    loop->LoadMagField("../data/elmag/VecB2.txt");

    loop->AddTask(new DriftTimeTask());
    loop->AddTask(new XZPlotTask());
    loop->AddTask(new XYPlotTask());
    loop->AddTask(new GraphResTask());
    loop->AddTask(new HistResTask());

    RecoPadsTask* t = new RecoPadsTask(); 
    loop->AddTask(t);
    loop->AddTask(new CircleAndRKFitTask(t));

    loop->RunSingleLoop();

    // loop->ResetTasks();
    // loop->LoadRK("../data/rk_tracks/rk_tracks.root");
    // loop->AddTask(new CircleFitEnergyTask());
    // loop->AddTask(new PlotSelectionTask());

    // gErrorIgnoreLevel = 6001;
    // loop->RunRKLoop();

    return 0;
}

int main()
{
    return reco_track();
}