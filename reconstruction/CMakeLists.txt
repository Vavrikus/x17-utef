# Check CMake version.
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(X17)

# ---Build executables--------------------------------------------------------------------
add_executable(reco_track reco_track.cpp)

target_link_libraries(reco_track Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)

target_include_directories(reco_track PRIVATE ../include)