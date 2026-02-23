# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(ionoel_pad_dist)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(ionoel_pad_dist ionoel_pad_dist.cpp)
target_link_libraries(ionoel_pad_dist Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
