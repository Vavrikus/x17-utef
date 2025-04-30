# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(thesis_plots)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(pad_layout pad_layout.cpp)
target_link_libraries(pad_layout Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
add_executable(test_track test_track.cpp)
target_link_libraries(test_track Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
