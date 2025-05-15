# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(mardia_test)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(mardia_test mardia_test.cpp)
target_link_libraries(mardia_test Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
