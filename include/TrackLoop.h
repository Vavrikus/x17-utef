#pragma once

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

#include "Track.h"
#include "VectorField.h"
#include "X17Utilities.h"

class TrackLoop;

/// @brief Interface for creating tasks for TrackLoop
class Task
{
    friend class TrackLoop;
protected:
    TrackLoop* loop = nullptr; // loop that will run this task
private:
    /// @brief To be run before the electron looping
    virtual void PreLoop() = 0;

    /// @brief To be run during the electron looping
    virtual void Loop() = 0;

    /// @brief To be run after the electron looping
    virtual void PostLoop() = 0;
};

/// @brief Task that includes plotting
class PlotTask : public Task
{
public:
    bool makeNewCanvas = true;    // should new canvas be made for this plot during postloop
    std::string canvasName = "";  // name of new canvas
    std::string canvasTitle = ""; // title of new canvas
    TCanvas* canvas = nullptr;    // pointer to canvas for plotting
};

/// @brief Struct for storing initial and final coordinates of ionization electron
struct IonElectron
{
    double x0,y0,z0; // Initial position [cm]
    double t0;       // Initial time [ns]
    double x1,y1,z1; // Final position [cm]
    double t1;       // Final time [ns]
};

/// @brief Class for iterating through tracks and performing different tasks
class TrackLoop
{   
private:
    Field<SensorData>* map;      // ionization electron map
    VectorField* magfield;       // magnetic field
    std::vector<Task*> tasklist; // contains all tasks to be run

    TTree* single_track = nullptr; // single track information (initial and final ionization electron positions)
    IonElectron curr_electron;     // variable for electron positions ([cm] or [ns])
    SensorData curr_reco;          // current reconstructed electron [cm]

    TTree* rk_tracks = nullptr; // information from Runge-Kutta tracks
    TrackRK* curr_track;         // current Runge-Kutta track

    /// @brief Reconstructs current electron coordinates
    void RecoElectron() {curr_reco = map->Invert(curr_electron.x1,curr_electron.y1,curr_electron.t1);}

    /// @brief To be run before the electron looping
    void PreLoop() {for (Task* t : tasklist) t->PreLoop();}

    /// @brief To be run during the electron looping
    void Loop()    {for (Task* t : tasklist) t->Loop();}

    /// @brief To be run after the electron looping
    void PostLoop() 
    {
        for (Task* t : tasklist) 
        {
            PlotTask* pt = dynamic_cast<PlotTask*>(t);
            if ((pt != nullptr) && pt->makeNewCanvas)
            {
                TCanvas* c = new TCanvas(pt->canvasName.c_str(),pt->canvasTitle.c_str());
                pt->canvas = c;
            }

            t->PostLoop();
        }
    }

public:
    /// @brief Destructor
    ~TrackLoop()
    {
        delete map;
        // for (Task* t : tasklist) delete t;
        if (single_track != nullptr) delete single_track;
    }

    /// @brief Adds task to TrackLoop task list
    /// @param t Task to be added
    void AddTask(Task* t)
    {
        t->loop = this;
        tasklist.push_back(t);
    }

    /// @brief  Getter for current electron information
    /// @return Reference to member curr_electron
    const IonElectron& GetCurrentElectron() {return curr_electron;}

    /// @brief  Getter for current reconstructed electron information
    /// @return Reference to member curr_reco
    const SensorData& GetCurrentReco() {return curr_reco;}

    /// @brief  Getter for current Runge-Kutta track
    /// @return Member curr_track
    const TrackRK* GetCurrentTrack() {return curr_track;}

    /// @brief  Getter for magnetic field
    /// @return Member magfield
    VectorField* GetMagField() {return magfield;}

    /// @brief  Getter for ionization electron map
    /// @return Member map
    Field<SensorData>* GetMap() {return map;}

    /// @brief  Getter for TTree with ionization electrons from single track
    /// @return Member single_track
    TTree* GetSingleTrack() {return single_track;}

    /// @brief Loads magnetic field from txt file (units = meters)
    /// @param path Path to the file
    void LoadMagField(const char* path, bool printInfo = false)
    {
        magfield = new VectorField(-0.2,0.2,-0.3,0.3,-0.3,0.3,0.005);
        magfield->LoadField(path);

        if (printInfo)
        {
            double minfield,maxfield,minangle,maxangle;
            X17::GetMinMaxField(*magfield,minfield,maxfield);
            X17::GetMinMaxFieldAngle(*magfield,minangle,maxangle);
            cout << "At least 0.0 cm from TPC walls: minimal magnetic field: " << minfield << " maximal: " << maxfield;
            cout << " minimal angle to electric field: " << minangle << " maximal: " << maxangle << "\n";
            X17::GetMinMaxField(*magfield,minfield,maxfield,0.5);
            X17::GetMinMaxFieldAngle(*magfield,minangle,maxangle,0.5);
            cout << "At least 0.5 cm from TPC walls: Minimal magnetic field: " << minfield << " maximal: " << maxfield;
            cout << " minimal angle to electric field: " << minangle << " maximal: " << maxangle << "\n";
            X17::GetMinMaxField(*magfield,minfield,maxfield,1);
            X17::GetMinMaxFieldAngle(*magfield,minangle,maxangle,1);
            cout << "At least 1.0 cm from TPC walls: Minimal magnetic field: " << minfield << " maximal: " << maxfield;
            cout << " minimal angle to electric field: " << minangle << " maximal: " << maxangle << "\n";
        }
    }

    /// @brief Loads the ionization electron map from a ROOT file
    /// @param path Path to the file
    void LoadMap(const char* path)
    {
        TFile* map_input = new TFile(path);
        map = (Field<SensorData>*)map_input->Get("map");
        // delete map_input;
    }

    /// @brief Loads file with Runge-Kutta simulated tracks
    /// @param path Path to file
    void LoadRK(const char* path)
    {
        TFile* input = new TFile(path);
        rk_tracks = (TTree*)input->Get("rk_tracks");
        // delete input;
    }

    /// @brief Loads single track simulation ROOT file
    /// @param path Path to file
    void LoadSingle(const char* path)
    {
        TFile* input = new TFile(path);
        single_track = (TTree*)input->Get("electrons");
        // delete input;
    }

    /// @brief Removes all tasks from task list
    void ResetTasks() {tasklist.clear();}

    /// @brief Runs all tasks for single track
    void RunSingleLoop()
    {
        PreLoop();

        // setting variables from single_track TTree
        single_track->SetBranchAddress("x0",&curr_electron.x0);
        single_track->SetBranchAddress("y0",&curr_electron.y0);
        single_track->SetBranchAddress("z0",&curr_electron.z0);
        single_track->SetBranchAddress("t0",&curr_electron.t0);
        single_track->SetBranchAddress("x1",&curr_electron.x1);
        single_track->SetBranchAddress("y1",&curr_electron.y1);
        single_track->SetBranchAddress("z1",&curr_electron.z1);
        single_track->SetBranchAddress("t1",&curr_electron.t1);
        
        // looping through all electrons
        int n_electrons = 0;
        for (int i = 0; i < single_track->GetEntries(); ++i)
        {
            single_track->GetEntry(i);
            if (X17::IsInSector(curr_electron.x1,curr_electron.y1,0) && X17::IsInSector(curr_electron.x0,curr_electron.y0,curr_electron.z0,-0.01)) 
            {
                n_electrons++;
                RecoElectron();
                Loop();
            }
        }
        cout << "\nNumber of electrons in the TPC region: " << n_electrons << "\n";

        PostLoop();
    }

    /// @brief Runs all tasks for Runge-Kutta simulated tracks
    void RunRKLoop()
    {
        PreLoop();

        rk_tracks->SetBranchAddress("track",&curr_track);
        int n_tracks = rk_tracks->GetEntries();

        for (int i = 0; i < n_tracks; i++)
        {
            rk_tracks->GetEntry(i);
            if((100*i)%n_tracks == 0) std::cout << "Progress: " << 100*i/n_tracks << " \%\n";
            std::cout << "Track " << i+1 << " out of " << n_tracks << ".\n";
            Loop();
            // if (i == 1999) break;
        }

        PostLoop();        
    }
};

#include "Tasks.h"