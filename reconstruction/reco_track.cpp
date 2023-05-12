// ROOT dependencies
#include "TError.h"
#include "TFile.h"
#include "TTree.h"

// X17 dependencies
#include "TrackLoop.h"

#include "RecoTasks.h"

using namespace X17;

int reco_track()
{
    // Loading the magnetic field data.
    X17::Field<X17::Vector>* magfield = X17::LoadField("../data/elmag/VecB2.txt",{-20,-30,-30},{20,30,30},0.5);

    // Loading the ionization electron drift map.
    TFile* map_input = new TFile("../data/ion_map/map.root");
    X17::Field<X17::MapPoint>* map = (X17::Field<X17::MapPoint>*)map_input->Get("map");

    // Loading file with one microscopic track.
    TFile* input = new TFile("../data/single_track/new_9010/electrons.root");
    TTree* single_track = (TTree*)input->Get("electrons");

    // Loading file with Runge-Kutta tracks.
    TFile* input2 = new TFile("../data/rk_tracks/rk_tracks.root");
    TTree* rk_tracks = (TTree*)input2->Get("rk_tracks");


    // TrackLoop for single microscopic track.
    TrackLoop* single_loop = new TrackLoop(map,magfield);
    
    single_loop->AddTask(new DriftTimeTask());
    single_loop->AddTask(new XZPlotTask());
    single_loop->AddTask(new XYPlotTask());
    single_loop->AddTask(new GraphResTask());
    single_loop->AddTask(new HistResTask());

    RecoPadsTask* t = new RecoPadsTask();
    single_loop->AddTask(t);
    single_loop->AddTask(new CircleAndRKFitTask(t));


    // TrackLoop for Runge-Kutta simulated tracks.
    TrackLoop* rk_loop = new TrackLoop(map,magfield);
    rk_loop->AddTask(new CircleFitEnergyTask());
    // rk_loop->AddTask(new PlotSelectionTask());

    gErrorIgnoreLevel = 6001;

    // Processing.
    // single_loop->ProcessSingle(single_track);
    rk_loop->ProcessRK(rk_tracks);

    return 0;
}

int main(int argc, char const *argv[])
{
    return reco_track();
}