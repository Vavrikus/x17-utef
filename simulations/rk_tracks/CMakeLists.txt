# Check CMake version.
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(X17)

# ---Build executables--------------------------------------------------------------------
add_executable(rk_tracks rk_tracks.cpp)

target_link_libraries(rk_tracks Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_include_directories(rk_tracks PRIVATE ../../include)