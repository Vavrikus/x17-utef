# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(Electron_positron)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(make_track make_track.cpp)
target_link_libraries(make_track Garfield::Garfield)
