# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(bfield)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(bfield bfield.C)
target_link_libraries(bfield Garfield::Garfield)
