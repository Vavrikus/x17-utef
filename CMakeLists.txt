# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(X17)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

add_compile_options(-O3) # enable maximum optimization

# ---Root---------------------------------------------------------------------------------
message(STATUS "Looking for ROOT")
find_package(ROOT REQUIRED)
message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_BINDIR}")

set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS}")

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

find_package(GSL REQUIRED)

# ---Add subdirectory containing source code----------------------------------------------
include_directories(include)
file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp")
add_library(X17_source_lib STATIC ${SOURCES})

# ---Build executables--------------------------------------------------------------------
add_executable(reco_track reconstruction/reco_track.cpp)
add_executable(ion_electrons simulations/ion_map/ion_electrons.cpp)
add_executable(make_map simulations/ion_map/make_map.cpp)
add_executable(rk_tracks simulations/rk_tracks/rk_tracks.cpp)
add_executable(make_track simulations/single_track/make_track.cpp)
add_executable(gas_table simulations/single_track/gas_table.cpp)

target_link_libraries(reco_track Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_link_libraries(ion_electrons Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_link_libraries(make_map Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_link_libraries(rk_tracks Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_link_libraries(make_track Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)
target_link_libraries(gas_table Garfield::Garfield ${ROOT_LIBRARIES} X17_source_lib)

target_include_directories(reco_track PRIVATE include)
target_include_directories(ion_electrons PRIVATE include)
target_include_directories(make_map PRIVATE include)
target_include_directories(rk_tracks PRIVATE include)
target_include_directories(make_track PRIVATE include)
target_include_directories(gas_table PRIVATE include)