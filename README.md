# X17 TPC simulation and reconstruction
C++ code in this project uses ROOT and Garfield++ to simulate and reconstruct electron and positron tracks inside TPC. The workflow is split between several scripts.

## Installation
The user needs to have ROOT installed (built on 6.26/06), for running simulations Garfield++ has to be installed. Currently used simulation results are provided. ROOT scripts can be run using the ROOT command:

```bash
root name_of_script.cpp
```
Simulation scripts have to be compiled before running:

```bash
cd /script_folder/build
cmake ..
make
./script_name
```

### Map simulation
The map of the ionization electron drift can be simulated on MetaCentrum:
1. Files need to be synchronized and compiled (can be achieved by adjusting *sync.sh* and *remote_compile.sh*).
2. Directory *data* should be cleared so no files get overwritten.
3. Simulation is started by executing *ion_multiple.sh*:

    ```bash
    ./ion_multiple.sh num_of_threads step_size_cm
    ```
    Make sure wall time and memory inside *ion_single.sh* are sufficient.
4. After job finishes, output files ion(id).root and ion(id).out are sent into the data directory. You can fetch them by adjusting *fetch_data.sh*.
5. Map can be compiled using script *make_map.cpp* (ion.root files have to be in the data directory).


## Usage
The project uses the following folder structure (headers and less important files ommited):
- **The project folder**
    - *GEMReadout_CenteredConector.kicad_pcb* -- KiCad file with the readout PCB
    - *tpc_layout.ggb* -- GeoGebra file with a top-down view of the TPCs
    - *sync.sh* -- Syncing with MetaCentrum (needs adjustments)
    - **.vscode**
        - *root-on-vscode.code-workspace* -- useful setup for using ROOT and Garfield++ in Visual Studio Code with Intellisense and debuging tool (some paths may need adjusting)
    - **Magnetic_field_plot**
        - *bfield.C* -- old script for displaying the magnetic field using Garfield++
    - **mag_data**
        - Contains files *VecB.txt* and *VecE.txt* with magnetic field simulation (electric field is homogeneous) and files *VecB2.txt* and *VecE2.txt* containing the same data with adjusted coordinates.
    - **presentations**
        - Contains presentation from 07-03-2023 seminar (tex and pdf files)
    - **electron_positron_tracks**
        - *circle_lines.ggb* -- GeoGebra file with circle fit geometry
        - *fetch_data.sh*, *ion_multi.sh*, *ion_single.sh*, *remote_compile.sh* -- scripts for handling the map simulation on MetaCentrum (see above)
        - *gas_table.cpp* -- script for generating the gas table for AvalancheMC simulation (not necessary for microscopic tracking)
        - *ion_electrons.cpp* -- script for ionization electrons simulation (creates input for *make_map.cpp*)
        - *make_map.cpp* -- script for creating the map of ionization electrons drift
        - *make_track.cpp* -- script for single track simulation
        - *reco_track.cpp* -- script for track reconstruction