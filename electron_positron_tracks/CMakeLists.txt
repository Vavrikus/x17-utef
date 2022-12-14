# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(Electron_positron)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Root------------------------------------------------------------------------
find_package(ROOT REQUIRED)

set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS}")

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

find_package(GSL REQUIRED)

# ---Build executable------------------------------------------------------------
add_executable(make_track make_track.cpp)
add_executable(ion_electrons ion_electrons.cpp)
add_executable(gas_table gas_table.cpp)

target_link_libraries(make_track Garfield::Garfield ${ROOT_LIBRARIES})
target_link_libraries(ion_electrons Garfield::Garfield ${ROOT_LIBRARIES})
target_link_libraries(gas_table Garfield::Garfield ${ROOT_LIBRARIES})