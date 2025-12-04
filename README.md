# X17 TPC simulation and reconstruction
C++ code in this project uses ROOT and Garfield++ to simulate and reconstruct electron and positron tracks inside the atypic TPC of the X17 project. The workflow is split between several scripts.

## Installation
The user needs to have ROOT installed (built on 6.26/06), for running simulations Garfield++ has to be installed. Currently used simulation results are provided.
Scripts have to be compiled before running:

```bash
cd /path/to/project/folder
mkdir build && cd build
cmake ..
make
./script_name
```

Alternatively, the script *rebuild_debug.sh* or *rebuild_release.sh* can be used. All of the binaries will be saved to the newly made build folder. This might lead to errors in the scripts whenever relative paths are used.

### Map simulation
The map of the ionization electron drift can be simulated on MetaCentrum:
1. Files need to be synchronized and compiled (can be achieved by adjusting *sync.sh* and *remote_compile.sh*).
2. Directory *data/ion_map/new_sample* should be cleared so no files get overwritten.
3. Simulation is started by executing *ion_multiple.sh*:

    ```bash
    ./ion_multiple.sh num_of_threads step_size_cm num_of_iterations_per_electron
    ```
    Make sure wall time and memory inside *ion_single.sh* are sufficient.
4. After the job finishes, output files ion(id).root and ion(id).out are sent into the data directory. You can fetch them by adjusting *fetch_data.sh*.
5. Map can be compiled using script *make_map.cpp* (ion.root files have to be in the data directory).

### Microscopic track simulation
Microscopic tracks can be simulated on MetaCentrum:
1. Files need to be synchronized and compiled (can be achieved by adjusting *sync.sh* and *remote_compile.sh*).
2. Directory *data/micro_tracks/new_tracks* should be cleared so no files get overwritten.
3. Simulation is started by executing *track_multi.sh*:

    ```bash
    ./track_multi.sh num_of_threads num_of_iterations_per_track num_of_angle_bins num_of_energy_bins
    ```
    If only *num_of_threads* parameter is provided the track parameters will be random instead of grid-like. Make sure wall time and memory inside *track_single.sh* are sufficient.
4. After the job finishes, output files track_full(id).root, track_small(id).root and tracks(id).out are sent into the data directory. You can fetch them by adjusting *fetch_data.sh*.
5. The tracks can be reconstructed with the reco_track script.

## Usage
The project uses the following folder structure:
- **The project folder**
    - **.vscode**
        - *root-on-vscode.code-workspace* and *launch.json* -- useful setup for using ROOT and Garfield++ in Visual Studio Code with Intellisense and debuging tool (include paths need to be adjusted)
    - **tests/magnetic_field_plot**
        - *bfield.C* -- old script for displaying the magnetic field using Garfield++
    - **data**
        - **elmag**
            - Contains files *VecB.txt* and *VecE.txt* with magnetic field simulation (electric field is homogeneous) and files *VecB2.txt* and *VecE2.txt* containing the same data with adjusted coordinates.
        - **ion_map**
            - Contains files with ionization electron data and compiled map in the file *map.root*
            - *sample_1.0* contains the old map with 90/10 gas composition
            - *sample_2.0* contains the new map with 70/30 gas composition
        - **micro_tracks**
            - Contains files with microscopic track simulations (many tracks)
            - *plot_drift.C* -- ROOT script for driftline plots used in the RD51 presentation
            - Currently used tracks are in *grid_01* and *grid_02*
        - **single_track**
            - Contains files with single simulated track for different settings
            - *original*, *original_v2*, *originalMC*, *track1* are tracks with 8 MeV momentum and use old coordinates
            - *new_9010*, *new_7030* are tracks with 8 MeV momentum and use new coordinates; the latter is the only track in this group with non-0.1 eV starting ionization electron energies
            - *newest_7030* is a track with 10 MeV (?) momentum and uses new coordinates and struct MicroPoint
    - **figures**
        - Contains the newest figures (post-thesis, to be used in the article, etc.)
    - **include**
        - Contains headers for the project with doxygen comments (and .inl files with definitions of templated functions)
    - **presentations**
        - Contains presentation from meetings and seminars (tex and pdf files)
    - **reconstruction**
        - *map_test.cpp* -- script for plots of the map reconstruction residues of microscopic tracks
        - *reco_mtracks.cpp* -- script for reconstruction of microscopic tracks
        - *reco_plots.cpp* -- script for plots of the reconstruction (comparison of simulated and reconstructed energy)
        - *reco_test.cpp* -- script for reconstruction testing (single track plots, Runge-Kutta tracks circle fit)
    - **remote**
        - Useful scripts for MetaCentrum (uses absolute paths, needs to be adjusted)
        - *sync.sh* -- Syncing with MetaCentrum
        - *clear_data.sh* -- Clearing the data directories on MetaCentrum
        - *fetch_data.sh* -- Syncing data simulated on MetaCentrum
        - *remote_compile.sh* -- Compiling scripts on MetaCentrum
    - **schematics**
        - *circle_fit2d.ggb* -- GeoGebra file with 2D circle fit geometry
        - *circle_lines.ggb* -- GeoGebra file with 3D circle fit geometry
        - *GEMReadout_CenteredConector.kicad_pcb* -- KiCad file with the readout PCB
        - *map_visualization.ggb* -- Simple visualization of the map and its inverse
        - *phi_systematic_error.ggb* -- Visualization of the phi of the maxima and minima of systematic error in reconstructed energy (pads shown)
        - *tpc_layout.ggb* -- GeoGebra file with a top-down view of the TPCs (area of the first map simulation and of the magnetic field simulation shown)
        - *tpc_micro_simulation* -- Visualization of the 3D layout of the microscopic track simulation
    - **simulations**
        - **ion_map**
            - *ion_electrons.cpp* -- script for ionization electrons simulation (creates input for *make_map.cpp*)
            - *make_map.cpp* -- script for creating the map of ionization electron drift
            - *ion_multi.sh*, *ion_single.sh* -- scripts handling the map simulation on MetaCentrum (see above)
        - **rk_tracks**
            - *rk_tracks.cpp* -- script for quick Runge-Kutta track simulation, generates rk_tracks.root file (currently 100,000 tracks)
        - **single_track**
            - *gas_table.cpp* -- script for generating the gas table for AvalancheMC simulation (not necessary for microscopic tracking)
            - *make_track.cpp* -- script for single track simulation
        - **micro_tracks**
            - *micro_tracks.cpp* -- script for random or grid-like simulation of multiple microscopic tracks
            - *track_multi.sh*, *track_single.sh* -- scripts handling the simulation on MetaCentrum
    - **source**
        - Contains all .cpp files common to the reconstruction and simulation scripts.
    - **thesis**
        - Contains the latex and pdf files of my bachelor thesis