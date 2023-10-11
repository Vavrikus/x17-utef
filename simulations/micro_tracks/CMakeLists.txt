# Check CMake version.
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(X17)

# find_package(GSL REQUIRED)
include_directories(/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/include)
link_directories(/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/lib)

# ---Build executables--------------------------------------------------------------------
add_executable(micro_tracks micro_tracks.cpp)

target_link_libraries(micro_tracks Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib ROOT::RIO)
target_include_directories(micro_tracks PRIVATE ../../include)