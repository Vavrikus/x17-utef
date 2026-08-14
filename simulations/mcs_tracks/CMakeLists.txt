# Check CMake version.
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(X17)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# find_package(GSL REQUIRED)
include_directories(/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/include)
link_directories(/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/lib)

# ---Build executables--------------------------------------------------------------------
add_executable(mcs_tracks mcs_tracks.cpp)

target_link_libraries(mcs_tracks Garfield::Garfield ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} X17_source ROOT::RIO)
target_include_directories(mcs_tracks PRIVATE ../../include)