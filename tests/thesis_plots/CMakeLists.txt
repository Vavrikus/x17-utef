# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(thesis_plots)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(thesis_plots plots.cpp)
target_link_libraries(thesis_plots Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
