# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(plot_drift)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(plot_drift plot_drift.cpp)
target_link_libraries(plot_drift Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
